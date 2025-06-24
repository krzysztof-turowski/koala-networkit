#pragma once

#include <boost/functional/hash.hpp>
#include <networkit/graph/Graph.hpp>
#include <vector>

using pairDistance_t =
    std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
class TopologyHeap {
   public:
    std::pair<int, NetworKit::node> top();

    void closeNode(NetworKit::node);

    bool empty();

    void printStorage();

    TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions, pairDistance_t& distances,
                 int source, const std::unordered_set<NetworKit::node>& extraBoundryNodes);

   private:
    void initializeStorage();
    void batchUpdate(int batch, NetworKit::node node, int currentValue);
    void fixIndex(int i);
    void initialDistances();

    static const int INF = std::numeric_limits<int>::max();  // TODO rewrite logic to handle double

    NetworKit::node source;
    NetworKit::node target;

    NetworKit::Graph& graph;
    nodeSubsets_t& regions;
    pairDistance_t& distances;
    std::vector<int> nodeStorageIdx;
    std::vector<int> isClosed;
    const std::unordered_set<NetworKit::node> extraBoundryNodes;
    std::vector<std::vector<int>> nodeRegions;
    std::vector<std::vector<int>> batches;
    std::vector<std::vector<int>> regionBatches;
    std::vector<std::vector<int>> sortedBoudryNodeRegions;

    int size;
    int storageSize;
    std::vector<std::pair<int, int>> storage;
};

} /* namespace Koala */
