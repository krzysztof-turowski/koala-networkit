#pragma once

#include <unordered_map>

#include "flow/MaximumFlow.hpp"

namespace Koala {

/**
 * @ingroup flow
 * Malhotra-Kumar-Maheshwari maximum-flow algorithm.
 *
 * The algorithm constructs a layered residual graph and computes each blocking flow by
 * repeatedly selecting a minimum-potential vertex and pushing flow forward and backward.
 *
 * @see V. M. Malhotra, M. Pramodh Kumar, and S. N. Maheshwari,
 *      "An O(|V|^3) Algorithm for Finding Maximum Flows in Networks", 1978.
 */
class MalhotraKumarMaheshwariFlow final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;

    /**
     * Compute a maximum flow with layered residual graphs and blocking flows.
     */
    void run();

 private:
    std::unordered_map<NetworKit::Edge, int> flow;
    NetworKit::Graph graph_stage;
    std::unordered_map<NetworKit::node, int> level;
    std::unordered_map<NetworKit::node, int> in_potential, out_potential;
    std::unordered_map<NetworKit::Edge, int> capacity;

    NetworKit::Edge reverse(const NetworKit::Edge &);

    bool build_level_graph();
    void compute_potential();
    void push_forward(NetworKit::node, NetworKit::edgeweight);
    void push_backward(NetworKit::node, NetworKit::edgeweight);
    void delete_node(NetworKit::node);
    void initialize();
};

}  // namespace Koala
