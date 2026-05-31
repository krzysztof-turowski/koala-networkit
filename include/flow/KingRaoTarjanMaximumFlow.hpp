#pragma once

#include "flow/PushRelabelMaximumFlow.hpp"
#include "flow/maximum_flow/KrtEdgeDesignator.hpp"

namespace Koala {

/**
 * @ingroup flow
 * The class for the King-Rao-Tarjan maximum flow algorithm
 */
class KingRaoTarjanMaximumFlow final : public PushRelabelMaximumFlow {
 public:
    KingRaoTarjanMaximumFlow(
        NetworKit::Graph&, NetworKit::node, NetworKit::node,
        KRTEdgeDesignator::Parameters = {});
    void run();

 private:
    KRTEdgeDesignator edge_designator;
    KRTEdgeDesignator::Parameters edge_designator_parameters;

    NetworKit::node get_active_vertex() override;
    NetworKit::node get_admissible_residual_edge(NetworKit::node) override;
    void on_relabel(NetworKit::node, int) override;
};

}  /* namespace Koala */
