#pragma once

#include "flow/PushRelabelMaximumFlow.hpp"
#include "flow/maximum_flow/KrtEdgeDesignator.hpp"

namespace Koala {

/**
 * @ingroup flow
 * King-Rao-Tarjan deterministic maximum-flow algorithm.
 *
 * This implementation combines the shared push-relabel engine with an edge designator that
 * supplies residual-edge candidates and responds to relabeling and failed candidates.
 *
 * @see https://doi.org/10.1006/jagm.1994.1044
 */
class KingRaoTarjanMaximumFlow final : public PushRelabelMaximumFlow {
 public:
    /**
     * Construct a King-Rao-Tarjan maximum-flow instance.
     *
     * @param graph Input graph.
     * @param source Source vertex.
     * @param target Sink vertex.
     * @param edge_designator_parameters Optional designator tuning parameters.
     */
    KingRaoTarjanMaximumFlow(
        NetworKit::Graph&, NetworKit::node, NetworKit::node,
        KRTEdgeDesignator::Parameters = {});

    /**
     * Initialize the designator and compute a maximum flow.
     */
    void run();

 private:
    KRTEdgeDesignator edge_designator;
    KRTEdgeDesignator::Parameters edge_designator_parameters;

    NetworKit::node get_active_vertex() override;
    NetworKit::node get_admissible_residual_edge(NetworKit::node) override;
    void on_relabel(NetworKit::node, int) override;
};

}  // namespace Koala
