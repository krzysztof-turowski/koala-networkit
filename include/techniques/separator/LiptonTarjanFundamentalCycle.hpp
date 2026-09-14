#pragma once

#include <optional>
#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

#include "shortest_path/planar/SuitableRDivision.hpp"

namespace Koala {
using cycle_t = std::vector<NetworKit::node>;
using index_map_t = std::vector<std::unordered_map<NetworKit::node, int>>;

struct CostComputation {
    double insideCost;
    double outsideCost;
    bool insideIsClockwiseArc;
};

class LiptonTarjanFundamentalCycle : public NetworKit::Algorithm {
   public:
    void run() override;
    std::optional<cycle_t> getFundamentalCycle();
    LiptonTarjanFundamentalCycle(NetworKit::Graph& G,
                                 planar_embedding_t& embedding,
                                 std::vector<NetworKit::node>& parent,
                                 std::vector<double>& costsOfSubtree,
                                 std::vector<double>& vertexCost,
                                 NetworKit::node root);

   private:
    void find_indexing_map();
    std::optional<cycle_t> build_fundamental_cycle(
        NetworKit::node v1, NetworKit::node w1,
        std::vector<NetworKit::node> parent);
    CostComputation compute_sides_initial_cost(const cycle_t& cycle);
    cycle_t shrink_fundamental_cycle(const cycle_t& initial_cycle,
                                     const CostComputation& sides_cost,
                                     const index_map_t& idx_of,
                                     NetworKit::node root, NetworKit::node v1,
                                     NetworKit::node w1);

    std::tuple<std::vector<NetworKit::node>, std::vector<bool>, NetworKit::node>
    compute_path_to_cycle(NetworKit::node apex, NetworKit::node curVi,
                          NetworKit::node curWi, std::vector<bool>& isOnCycle);
    NetworKit::node getApex(NetworKit::node v, NetworKit::node u, int dir);

    // helpers
    static NetworKit::node find_node_in_cycle_at_ith_position(
        const cycle_t& cycle, int i);
    static int wrapIndex(int i, int size);

    NetworKit::Graph& G;
    planar_embedding_t& embedding;
    std::vector<NetworKit::node>& parent;
    std::vector<double>& costs_of_subtree;
    std::vector<double>& vertex_cost;
    NetworKit::node root;
    std::optional<cycle_t> fundamental_cycle;
    index_map_t idx_of;
};
}  // namespace Koala
