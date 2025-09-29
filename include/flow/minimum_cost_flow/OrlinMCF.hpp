#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Attributes.hpp>

namespace Koala {

class OrlinMCF final : public MinimumCostFlow {
 public:
    // flow
    OrlinMCF(NetworKit::Graph const& g, edge_map<long long> const& costs,
        NetworKit::node s, NetworKit::node t, long long f) : MinimumCostFlow(g, costs, s, t, f) {
        init_flow();
    }
    // flow
    OrlinMCF(NetworKit::Graph const& g, edge_map<long long> const& costs,
            node_map<long long> const& b) : MinimumCostFlow(g, costs, b) {
        init_flow();
    }
    // circulation
    OrlinMCF(NetworKit::Graph const& g, edge_map<std::pair<long long,long long>> const& bounds,
            edge_map<long long> const& costs) : MinimumCostFlow(g, bounds, costs) {
        init_circulation();
    }

 private:
    void make_costs_nonnegative();

    const double ALPHA = 0.9;
    long long delta = 0;

    void run_impl() override;
    void init_flow() override;
    void init_circulation() override;
    long long cp(NetworKit::node, NetworKit::node);

    NetworKit::Graph generate_graph_for_sp(std::map<edge, std::pair<int, node>>&);
    bool is_imbalanced();
    void initialize();

    std::stack<std::pair<NetworKit::node, NetworKit::node>> contractions;
    void contract_nodes(NetworKit::node, NetworKit::node);
    void delta_scaling_phase();
    void contraction_phase();
    void uncontract_and_compute_flow();
    void scaling_augmentation(NetworKit::node, NetworKit::node);
    void fill_distances_uncapacitated(std::vector<double>&);
    NetworKit::Graph original_graph;
};

} /* namespace Koala */

 