#pragma once

#include <cstdint>
#include <unordered_map>
#include <utility>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <flow/minimum_cost_flow/MCFlowNetwork.hpp>

namespace Koala {

class MinimumCostFlow : public NetworKit::Algorithm {
 public:
    explicit MinimumCostFlow(MCFlowNetwork network) : network(network) {}
    
    void run() override {
        hasRun = false;
        run_impl();
        hasRun = true;
    }

    int64_t getMinCost() const {
        return min_cost;
    }

    std::unordered_map<NetworKit::Edge, int64_t> getMinCostFlow() const {
        return computed_flow;
    }

    int64_t getFlow(NetworKit::Edge const& edge) {
        if (computed_flow.find(edge) != computed_flow.end()) {
            return computed_flow.at(edge);
        }
        return 0;
    }

 protected:
    virtual void run_impl() = 0;
    MCFlowNetwork network;
    int64_t min_cost{0};
    std::unordered_map<NetworKit::Edge, int64_t> computed_flow;
};

} /* namespace Koala */
