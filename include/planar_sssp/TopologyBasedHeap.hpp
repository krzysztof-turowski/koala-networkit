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

        TopologyHeap(NetworKit::Graph& graph, nodeSubsets_t& regions, std::vector<pairDistance_t>& distances, int source);

    private:
        NetworKit::node source;
        NetworKit::node target;
        NetworKit::count distanceToTarget = 0;
        NetworKit::Graph& graph;
        nodeSubsets_t regions;
        std::vector<pairDistance_t> distances;
        int size;
        std::vector <std::pair<NetworKit::count, NetworKit::node>> storage;
    };

} /* namespace Koala */
