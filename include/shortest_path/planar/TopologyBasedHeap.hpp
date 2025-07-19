#pragma once

#include <networkit/graph/Graph.hpp>
#include <vector>

// Based on boost::hash_combine
struct pair_hash {
    size_t operator()(const std::pair<NetworKit::node, NetworKit::node>& pair) const {
        return pair.first ^
               (pair.second + 0x9e3779b97f4a7c15 + (pair.first << 12) + (pair.first >> 4));
    };
};

using pair_distance_t = std::unordered_map<std::pair<NetworKit::node, NetworKit::node>,
    NetworKit::edgeweight, pair_hash>;
using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
class TopologyHeap {
 public:
    std::pair<NetworKit::edgeweight, NetworKit::node> top();

    void close(NetworKit::node);

    bool empty();

    TopologyHeap(NetworKit::Graph& graph, node_subsets_t& regions, pair_distance_t& distances,
        NetworKit::node source, const std::unordered_set<NetworKit::node>& extra_boundary_nodes);

 private:
    void initialize_storage();
    void batch_update(
        NetworKit::index batch_id, NetworKit::node node, NetworKit::edgeweight current_value);
    void fix_index(NetworKit::index i);
    void initial_distances();

    NetworKit::node source;
    NetworKit::node target;

    NetworKit::Graph& graph;
    node_subsets_t& regions;
    pair_distance_t& distances;
    std::vector<NetworKit::index> node_storage_index;
    std::vector<int> is_closed;
    const std::unordered_set<NetworKit::node> extra_boundary_nodes;
    std::vector<std::vector<NetworKit::index>> node_regions;
    std::vector<std::vector<NetworKit::index>> batches;
    std::vector<std::vector<NetworKit::index>> region_batches;
    std::vector<std::vector<NetworKit::index>> sorted_boundary_node_regions;

    NetworKit::index size;
    NetworKit::index storage_size;
    std::vector<std::pair<NetworKit::edgeweight, NetworKit::node>> storage;
};

} /* namespace Koala */
