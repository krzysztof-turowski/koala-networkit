#pragma once

#include <boost/functional/hash.hpp>
#include <networkit/graph/Graph.hpp>
#include <vector>

using pair_distance_t =
    std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;
using node_subsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
class TopologyHeap {
 public:
    std::pair<int, NetworKit::node> top();

    void close(NetworKit::node);

    bool empty();

    TopologyHeap(NetworKit::Graph& graph, node_subsets_t& regions, pair_distance_t& distances,
        int source, const std::unordered_set<NetworKit::node>& extra_boundary_nodes);

 private:
    void initialize_storage();
    void batch_update(int batch, NetworKit::node node, int current_value);
    void fix_index(int i);
    void initial_distances();

    NetworKit::node source;
    NetworKit::node target;

    NetworKit::Graph& graph;
    node_subsets_t& regions;
    pair_distance_t& distances;
    std::vector<int> node_storage_index;
    std::vector<int> is_closed;
    const std::unordered_set<NetworKit::node> extra_boundary_nodes;
    std::vector<std::vector<int>> node_regions;
    std::vector<std::vector<int>> batches;
    std::vector<std::vector<int>> region_batches;
    std::vector<std::vector<int>> sorted_boundary_node_regions;

    int size;
    int storage_size;
    std::vector<std::pair<int, int>> storage;
};

} /* namespace Koala */
