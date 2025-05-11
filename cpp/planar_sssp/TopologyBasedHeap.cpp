#include "planar_sssp/TopologyBasedHeap.hpp"

#include <networkit/graph/Graph.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <vector>

namespace Koala {
    using pairDistance_t = std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;

    //root is stored in storage[1] to simplify traverse logic

    void TopologyHeap::fixIndex(int i) {
        int l = i * 2;
        int r = i * 2 + 1;

        if (storage[l].first < storage[r].first) {
            storage[i] = storage[l];
        }
        else {
            storage[i] = storage[r];
        }
    }

    void TopologyHeap::initialDistances() {
        //initial distances to boundry nodes from source

        if (nodeRegions[source].size() > 1) {//boundry node distances are already calculated
            storage[nodeStorageIdx[source]].first = 0;
            return;
        }

        int initialRegionId = nodeRegions[source][0];

        auto& region = regions[initialRegionId];
        auto subGraph = NetworKit::GraphTools::subgraphFromNodes(graph,
            std::unordered_set<NetworKit::node> {region.begin(), region.end()});
        std::vector<NetworKit::node> boundry;

        for (int batchId : regionBatches[initialRegionId]) {
            for (auto b : batches[batchId]) {
                boundry.push_back(b);
            }
        }
        NetworKit::MultiTargetDijkstra dij(subGraph, source, boundry.begin(), boundry.end());
        dij.run();
        auto map = dij.getTargetIndexMap();
        auto distances = dij.getDistances();
        for (auto b : boundry) {
            storage[nodeStorageIdx[b]].first = distances[map[b]];
        }
    }

    void TopologyHeap::initializeStorage() {

        // initialize the heap to have boundry sets in clusters
        std::vector<int> lastbatch;
        for (int i = 0; i < sortedBoudryNodeRegions.size(); i++) {
            int nodeIndex = sortedBoudryNodeRegions[i].back();
            sortedBoudryNodeRegions[i].pop_back();

            storage[i + storageSize] = { INF, nodeIndex };
            nodeStorageIdx[nodeIndex] = i + storageSize;

            if (lastbatch == sortedBoudryNodeRegions[i]) {
                batches.back().push_back(nodeIndex);
            }
            else {
                lastbatch = sortedBoudryNodeRegions[i];
                batches.push_back({ i });
            }
        }

        initialDistances();

        // build the heap
        for (int i = storageSize; i >= 1; i--) {
            fixIndex(i);
        }
    }

    TopologyHeap::TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions, pairDistance_t& distances, int source) :
        graph(graph), regions(regions), distances(distances), source(source) {

        //get list per node of regions a node is a part of
        nodeRegions.assign(graph.numberOfNodes(), std::vector<int>());
        for (int i = 0; i < regions.size(); i++) {
            for (auto node : regions[i]) {
                nodeRegions[node].push_back(i);
            }
        }

        // filter boundry nodes and add nodeId to the end
        for (int i = 0; i < graph.numberOfNodes(); i++) {
            if (nodeRegions[i].size() > 1) {
                nodeRegions[i].push_back(i);
                sortedBoudryNodeRegions.push_back(std::move(nodeRegions[i]));
            }
        }

        size = sortedBoudryNodeRegions.size();
        std::sort(sortedBoudryNodeRegions.begin(), sortedBoudryNodeRegions.end()); //TODO change sort from std::sort to some linera equivalent

        storageSize = 1;
        while (storageSize < size) storageSize <<= 1;

        //double the storageSize. Actual storage and heap structure
        storage.assign(storageSize << 1, { INF, -1 });
        initializeStorage();
    }

    std::pair<NetworKit::count, NetworKit::node> TopologyHeap::top() {
        assert(storage[0].second != -1);
        return storage[0];
    }

    void TopologyHeap::batchUpdate(int batchId, NetworKit::node node) {
        auto batch = batches[batchId];
        int nodeId = nodeStorageIdx[node];
        int minId = INF;
        std::deque<int> indexToFix;
        for (int b : batch) {
            int sId = nodeStorageIdx[b];
            if (indexToFix.back() != sId >> 1) {
                indexToFix.emplace_back(sId >> 1);
            }
            if (b == node) {
                storage[sId].first = INF;
            }
            else {
                storage[sId].first = std::min(storage[sId].first, storage[nodeId].first + distances[{node, b}]);
            }
        }

        while (!indexToFix.empty()) {
            int id = indexToFix.front();
            indexToFix.pop_front();

            if (id == 0) break;

            fixIndex(id);

            if (indexToFix.back() != id >> 1) {
                indexToFix.emplace_back(id >> 1);
            }
        }
    }

    void TopologyHeap::closeNode(NetworKit::node node) {
        for (int regionId : nodeRegions[node]) {
            for (int batchId : regionBatches[regionId]) {
                batchUpdate(batchId, node);
            }
        }
        size--;
    }

    bool TopologyHeap::empty() {
        return size == 0;
    }

} // namespace Koala
