#pragma once

#include <map>
#include <unordered_map>

#include "flow/MaximumFlow.hpp"

namespace Koala {

class PushRelabelMaximumFlow : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;
    void run();

 protected:
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> capacity, flow;
    std::map<NetworKit::node, int> distance, excess;

    virtual NetworKit::node get_active_vertex() = 0;
    virtual NetworKit::node get_admissible_residual_edge(NetworKit::node) = 0;
    virtual void on_relabel(NetworKit::node, int) { }

 private:
    void initialize();
    void push(NetworKit::node, NetworKit::node);
    void relabel(NetworKit::node);
};

}  /* namespace Koala */
