#pragma once

#include "flow/PushRelabelMaximumFlow.hpp"

namespace Koala {

class GoldbergTarjanPushRelabelMaximumFlow final : public PushRelabelMaximumFlow {
 public:
    using PushRelabelMaximumFlow::PushRelabelMaximumFlow;

 private:
    NetworKit::node get_active_vertex() override;
    NetworKit::node get_admissible_residual_edge(NetworKit::node) override;
};

}  /* namespace Koala */
