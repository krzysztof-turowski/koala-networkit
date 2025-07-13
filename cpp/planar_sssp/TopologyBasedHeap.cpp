#include "planar_sssp/TopologyBasedHeap.hpp"

#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <vector>

namespace Koala {
static const int INF = std::numeric_limits<int>::max();
using pairDistance_t =
    std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;

// root is stored in storage[1] to simplify traverse logic
void TopologyHeap::fixIndex(int i) {
    int l = i * 2, r = i * 2 + 1;

    if (storage[l].first < storage[r].first) {
        storage[i] = storage[l];
    } else {
        storage[i] = storage[r];
    }
}

void TopologyHeap::initialDistances() {
    // initial distances to boundary nodes from source

    if (nodeRegions[source].size() > 1 ||
        (extraBoundaryNodes.find(source) !=
            extraBoundaryNodes.end())) {  // boundary node distances are already calculated
        storage[node_storage_index[source]].first = 0;
    }

    int initialRegionId = nodeRegions[source][0];

    auto& region = regions[initialRegionId];
    auto subGraph = NetworKit::GraphTools::subgraphFromNodes(
        graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
    std::vector<NetworKit::node> boundary;

    for (int batchId : regionBatches[initialRegionId]) {
        for (auto b : batches[batchId]) {
            boundary.push_back(b);
        }
    }
    NetworKit::MultiTargetDijkstra dijkstra(subGraph, source, boundary.begin(), boundary.end());
    dijkstra.run();
    auto map = dijkstra.getTargetIndexMap();
    auto distances = dijkstra.getDistances();
    for (auto b : boundary) {
        storage[node_storage_index[b]].first = distances[map[b]];
    }
}

void TopologyHeap::initializeStorage() {
    // initialize the heap to have boundary sets in clusters
    std::vector<int> lastbatch;
    regionBatches.resize(regions.size());
    for (int i = 0; i < sortedBoudryNodeRegions.size(); i++) {
        int nodeIndex = sortedBoudryNodeRegions[i].back();
        sortedBoudryNodeRegions[i].pop_back();

        storage[i + storageSize] = {INF, nodeIndex};
        node_storage_index[nodeIndex] = i + storageSize;

        if (lastbatch == sortedBoudryNodeRegions[i]) {
            batches.back().push_back(nodeIndex);
        } else {
            lastbatch = sortedBoudryNodeRegions[i];
            batches.push_back({nodeIndex});
            for (auto region : sortedBoudryNodeRegions[i]) {
                regionBatches[region].push_back(batches.size() - 1);
            }
        }
    }

    initialDistances();

    // build the heap
    for (int i = storageSize - 1; i >= 1; i--) {
        fixIndex(i);
    }
}

TopologyHeap::TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions,
    pairDistance_t& distances, int source,
    const std::unordered_set<NetworKit::node>& extraBoundaryNodes)
    : graph(graph),
      regions(regions),
      distances(distances),
      source(source),
      extraBoundaryNodes(extraBoundaryNodes) {
    // get list per node of regions a node is a part of
    nodeRegions.assign(graph.numberOfNodes(), std::vector<int>{});
    node_storage_index.resize(graph.numberOfNodes());
    isClosed.assign(graph.numberOfNodes(), 0);
    for (int i = 0; i < regions.size(); i++) {
        for (auto& node : regions[i]) {
            nodeRegions[node].push_back(i);
        }
    }

    // filter boundary nodes and add nodeId to the end
    for (int i = 0; i < graph.numberOfNodes(); i++) {
        if (nodeRegions[i].size() > 1) {
            nodeRegions[i].push_back(i);
            sortedBoudryNodeRegions.push_back(nodeRegions[i]);
            nodeRegions[i].pop_back();
        } else if (extraBoundaryNodes.find(i) != extraBoundaryNodes.end()) {
            sortedBoudryNodeRegions.push_back({nodeRegions[i][0], i});
        }
    }

    size = sortedBoudryNodeRegions.size();
    std::sort(sortedBoudryNodeRegions.begin(), sortedBoudryNodeRegions.end());
    // TODO(289Adam289): change sort from std::sort to some linear equivalent

    storageSize = 1;
    while (storageSize < size) storageSize <<= 1;

    // double the storageSize. Actual storage and heap structure
    storage.assign(storageSize << 1, {INF, -1});

    initializeStorage();
}

std::pair<int, NetworKit::node> TopologyHeap::top() {
    if (size > 0) {
        assert(storage[1].second != -1);
    }
    return storage[1];
}

void TopologyHeap::batchUpdate(int batchId, NetworKit::node node, int currentValue) {
    auto batch = batches[batchId];
    int nodeId = node_storage_index[node];
    int minId = INF;
    std::deque<int> indexToFix;
    for (int b : batch) {
        int sId = node_storage_index[b];
        if (indexToFix.empty() || indexToFix.back() != sId >> 1) {
            indexToFix.emplace_back(sId >> 1);
        }
        if (isClosed[b]) continue;
        if (b == node) {
            continue;
        } else {
            int distance = INF;

            if (distances.find({node, b}) != distances.end()) {
                distance = currentValue + distances[{node, b}];
            }

            storage[sId].first = std::min(storage[sId].first, distance);
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
    size--;
    assert(nodeRegions[node].size() > 1 ||
           (extraBoundaryNodes.find(source) != extraBoundaryNodes.end()));

    int currentValue = storage[node_storage_index[node]].first;
    storage[node_storage_index[node]].first = INF;
    isClosed[node] = true;
    std::unordered_set<int> batcherToUpdate;
    for (int regionId : nodeRegions[node]) {
        for (int batchId : regionBatches[regionId]) {
            batcherToUpdate.insert(batchId);
        }
    }
    for (auto bId : batcherToUpdate) {
        batchUpdate(bId, node, currentValue);
    }
}

bool TopologyHeap::empty() {
    return size <= 0;
}

}  // namespace Koala
