#pragma once

#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * Residual-flow utility used by ElectricalFlow.
 *
 * The class stores a primal flow matrix, exposes residual capacity bounds, and supports
 * augmenting and rounding the computed flow.
 */
class FlowNetwork {
 public:
    /**
     * Construct an initially empty flow on a capacitated graph.
     *
     * @param graph Input capacitated graph.
     */
    explicit FlowNetwork(const NetworKit::Graph &graph);

    /**
     * @return The current total flow value.
     */
    double size() const;

    /**
     * @return The lower residual-capacity bound for edge (u, v).
     */
    double lowerCapacity(NetworKit::node u, NetworKit::node v) const;

    /**
     * @return The upper residual-capacity bound for edge (u, v).
     */
    double upperCapacity(NetworKit::node u, NetworKit::node v) const;

    /**
     * Round a nearly integral feasible flow while preserving its value.
     */
    void roundFlow();

    /**
     * Attempt to route an additional value from s to t in the residual graph.
     *
     * @return Whether the complete requested value was routed.
     */
    bool pushValue(NetworKit::node s, NetworKit::node t, double f);

    const NetworKit::Graph &graph;
    std::vector<std::vector<double>> flow;
    const int N, M;
};

}  // namespace Koala
