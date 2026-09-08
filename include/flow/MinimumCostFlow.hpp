#pragma once

#include <map>
#include <utility>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <flow/minimum_cost_flow/MCFlowNetwork.hpp>

namespace Koala {

class MinimumCostFlow : public NetworKit::Algorithm {
 public:
    explicit MinimumCostFlow(MCFlowNetwork const& network) : network(network) {}
    void run() {
        hasRun = false;
        run_impl();
        hasRun = true;
    }

    virtual std::int64_t getFlow(NetworKit::Edge const&) = 0;
    virtual std::unordered_map<NetworKit::Edge, std::int64_t> getMinCostFlow() const = 0;

    std::int64_t getMinCost() const {
        return min_cost;
    }

 protected:
    virtual void run_impl() = 0;
    MCFlowNetwork network;
    std::int64_t min_cost{0};
};

} /* namespace Koala */
