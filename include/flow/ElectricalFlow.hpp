#pragma once

#include <vector>

#include <networkit/graph/Graph.hpp>

#include "flow/MaximumFlow.hpp"
#include "flow/electrical_flow/FlowNetwork.hpp"

namespace Koala {

/**
 * @ingroup flow
 * Electrical-flow maximum-flow algorithm following Madry's augmenting-flow framework.
 *
 * The implementation routes electrical flows in a residual network, applies correction
 * steps, and optionally rounds the final primal flow. Directed inputs are reduced to an
 * initialized undirected instance internally.
 *
 * @see https://arxiv.org/abs/1608.06016
 */
class ElectricalFlow final : public MaximumFlow {
 public:
    /**
     * Construct an electrical-flow maximum-flow instance.
     *
     * @param graph Input graph.
     * @param s Source vertex.
     * @param t Sink vertex.
     * @param round Whether to round the final flow to an integral flow.
     */
    ElectricalFlow(
        NetworKit::Graph graph, NetworKit::node s, NetworKit::node t, bool round = true);

    /**
     * Compute a maximum flow with electrical-flow correction steps.
     */
    void run() override;

    /**
     * @return The original input graph.
     */
    const NetworKit::Graph &getGraph() const { return original_graph; }

    /**
     * @return The computed flow matrix indexed by edge endpoints.
     */
    const std::vector<std::vector<double>> &getFlow() const { return flow; }

 private:
    bool route_flow();
    void initialize();
    bool is_feasible();
    bool augmentation_step();
    void fixing_step();

    NetworKit::Graph original_graph;
    NetworKit::Graph graph;
    const NetworKit::node s, t;
    int U;
    int initial_flow;
    bool directed, round;

    std::vector<std::vector<double>> flow;
    std::vector<double> demand;
    FlowNetwork primal;
    std::vector<double> dual;
    double progress;
    double target_flow;
};

}  // namespace Koala
