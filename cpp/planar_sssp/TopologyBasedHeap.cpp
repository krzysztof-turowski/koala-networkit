#include "planar_sssp/TopologyBasedHeap.hpp"

#include <networkit/graph/Graph.hpp>
#include <vector>

namespace Koala {
    using pairDistance_t = std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;

    TopologyHeap::TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions, std::vector<pairDistance_t>& distances, int source) :
        graph(graph), regions(regions), distances(distances), source(source) {

    }

    std::pair<NetworKit::count, NetworKit::node> TopologyHeap::top() {

    }

    void TopologyHeap::closeNode(NetworKit::node) {

    }

    bool TopologyHeap::empty() {
        return size == 0;
    }

} // namespace Koala
