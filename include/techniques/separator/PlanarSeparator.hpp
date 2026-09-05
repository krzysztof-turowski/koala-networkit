#pragma once

#include <tuple>
#include <vector>

#include <graph/PlanarGraphTools.hpp>
#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#include "shortest_path/planar/SuitableRDivision.hpp"
#include "techniques/separator/BalancedSeparator.hpp"

namespace Koala {

class PlanarSeparator : public BalancedSeparator {
   public:
    explicit PlanarSeparator(const NetworKit::Graph& graph);
    PlanarSeparator(const NetworKit::Graph& graph,
                    const std::vector<double>& vertexCost);
    void run() override;

   private:
    std::vector<double> vertex_cost;
    void clean_partitions();

    std::pair<std::vector<int>, std::vector<NetworKit::node>>
    perform_BFS_and_find_spanning_tree(const NetworKit::Graph& G,
                                       NetworKit::node startNode);
    static std::vector<NetworKit::count> find_number_of_vertices_at_level(
        const std::vector<int>& lvl);
    std::pair<NetworKit::Graph, NetworKit::node> find_contracted_subgraph(
        const NetworKit::Graph& G, planar_embedding_t& embedding,
        std::vector<int>& lvl, int l0, int l2);
    std::tuple<int, int, int> find_partition_levels(
        const NetworKit::Graph& G,
        std::vector<NetworKit::node>& verticesAtLevel, std::vector<int>& lvl);

    std::vector<double> find_subtree_costs(
        const NetworKit::Graph& H, const std::vector<NetworKit::node>& parentH,
        const std::vector<double>& vertex_cost, NetworKit::node root);
    bool are_connected_components_eligible_for_partition(
        std::vector<std::vector<NetworKit::node>>& components);
    void find_separator_from_components(
        std::vector<std::vector<NetworKit::node>>& components);
    void extract_separator_and_partition_from_cycle(
        const NetworKit::Graph& G, const std::vector<int>& lvl, int l1, int l2,
        std::vector<NetworKit::node>& finalCycle, NetworKit::node x);
    void extract_level_partition(const NetworKit::Graph& G,
                                 const std::vector<int>& lvl, int l1, int l2);
    double compute_middle_cost(const NetworKit::Graph& G,
                               const std::vector<int>& lvl, int l1, int l2);
};

}  // namespace Koala
