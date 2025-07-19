#pragma once

#include "PlanarSSSP.hpp"

using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

class HenzingerPlanarSSSP : public PlanarSSSP {
 private:
    class HenzingerPriorityQueue {
        std::set<std::pair<NetworKit::count, NetworKit::node>> set;
        std::unordered_map<NetworKit::node, NetworKit::count> id_map;

        // caching the minimal element to ensure constant time read
        std::pair<NetworKit::count, NetworKit::node> minimum_element;

     public:
        void push(NetworKit::node id, NetworKit::count key);
        void update(NetworKit::node id, NetworKit::count new_key);
        void deactivate(NetworKit::node id);
        NetworKit::count minimum_key();
        NetworKit::node minimum_item();
        bool empty();
    };

    NetworKit::Graph normal_graph;
    std::vector<NetworKit::edgeweight> d;
    HenzingerPriorityQueue main_Q;
    std::vector<HenzingerPriorityQueue> Q;
    std::vector<std::vector<NetworKit::node>> regions;
    std::vector<int> is_boundary;
    NetworKit::count number_of_regions;
    int r;

    void initialize_queues(node_subsets_t& division);
    void main_thrust();

 public:
    HenzingerPlanarSSSP(
        NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target = NetworKit::none)
        : PlanarSSSP(graph, source, target) {}

    void run();
};

} /* namespace Koala */
