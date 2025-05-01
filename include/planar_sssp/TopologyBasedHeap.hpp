#pragma once

#include <networkit/graph/Graph.hpp>
#include <boost/functional/hash.hpp>
#include <vector>

using pairDistance_t = std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
    class TopologyHeap {
    public:

        std::pair<NetworKit::count, NetworKit::node> top();

        void closeNode(NetworKit::node);

        bool empty();

        TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions, pairDistance_t& distances, int source);

    private:
        void initializeStorage();
        void batchUpdate(int batch, NetworKit::node node);
        void fixIndex(int i);
        void initialDistances();

        static const int INF = std::numeric_limits<int>::max(); //TODO rewrite logic to handle double

        NetworKit::node source;
        NetworKit::node target;

        NetworKit::Graph& graph;
        nodeSubsets_t& regions;
        pairDistance_t& distances;
        std::vector<int> nodeStorageIdx;
        std::vector<std::vector<int>> nodeRegions;
        std::vector<std::vector<int>> batches;
        std::vector<std::vector<int>> regionBatches;
        std::vector<std::vector<int>> sortedBoudryNodeRegions;

        int size;
        int storageSize;
        std::vector <std::pair<NetworKit::count, NetworKit::node>> storage;
    };

} /* namespace Koala */
