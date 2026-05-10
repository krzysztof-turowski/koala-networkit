#pragma once

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {
class PlanarSSSP : public NetworKit::Algorithm {
 public:
    virtual void run() = 0;

    const std::vector<NetworKit::edgeweight> getDistances();
    NetworKit::edgeweight distance(NetworKit::node t);

    explicit PlanarSSSP(
        NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target = NetworKit::none);

 protected:
    NetworKit::Graph &graph;
    NetworKit::node source, target;
    std::vector<NetworKit::edgeweight> distances;
};

class FredericksonPlanarSSSP : public PlanarSSSP {
 private:
    NetworKit::Graph normal_graph;
    int c = 6;
    NetworKit::count r1 = NetworKit::none;
    NetworKit::count r2 = NetworKit::none;

 public:
    FredericksonPlanarSSSP(NetworKit::Graph &graph, NetworKit::node source, NetworKit::node target)
        : PlanarSSSP(graph, source, target) {}

    void setDivisionParameters(NetworKit::count level_1, NetworKit::count level_2);

    void run();
};

using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

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
    NetworKit::count r = NetworKit::none;

    void initialize_queues(node_subsets_t& division);
    void main_thrust();

 public:
    HenzingerPlanarSSSP(
        NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target = NetworKit::none)
        : PlanarSSSP(graph, source, target) {}

    void setDivisionParameter(NetworKit::count r);

    void run();
};

} /* namespace Koala */
