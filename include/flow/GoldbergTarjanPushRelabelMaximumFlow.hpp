#pragma once

#include "flow/PushRelabelMaximumFlow.hpp"

namespace Koala {

/**
 * @ingroup flow
 * Goldberg-Tarjan push-relabel maximum-flow algorithm.
 *
 * This implementation uses the common push-relabel engine with direct scans for active
 * vertices and admissible residual edges.
 *
 * @see https://doi.org/10.1145/48014.61051
 */
class GoldbergTarjanPushRelabelMaximumFlow final : public PushRelabelMaximumFlow {
 public:
    using PushRelabelMaximumFlow::PushRelabelMaximumFlow;

 private:
    NetworKit::node get_active_vertex() override;
    NetworKit::node get_admissible_residual_edge(NetworKit::node) override;
};

}  // namespace Koala
