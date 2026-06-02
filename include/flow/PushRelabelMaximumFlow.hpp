#pragma once

#include <map>
#include <unordered_map>

#include "flow/MaximumFlow.hpp"

namespace Koala {

/**
 * @ingroup flow
 * Common push-relabel maximum-flow engine.
 *
 * Push-relabel algorithms maintain a preflow and repeatedly push excess along admissible
 * residual edges or relabel active vertices. Concrete algorithms specialize how active
 * vertices and admissible residual edges are selected.
 *
 * @see https://doi.org/10.1145/48014.61051
 */
class PushRelabelMaximumFlow : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;

    /**
     * Compute a maximum flow using the subclass selection policy.
     */
    void run();

 protected:
    std::unordered_map<NetworKit::Edge, int> capacity, flow;
    std::map<NetworKit::node, int> distance, excess;

    virtual NetworKit::node get_active_vertex() = 0;
    virtual NetworKit::node get_admissible_residual_edge(NetworKit::node) = 0;
    virtual void on_relabel(NetworKit::node, int) { }

 private:
    void initialize();
    void push(NetworKit::node, NetworKit::node);
    void relabel(NetworKit::node);
};

}  // namespace Koala
