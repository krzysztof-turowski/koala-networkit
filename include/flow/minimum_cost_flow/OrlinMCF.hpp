#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Attributes.hpp>

namespace Koala {

class OrlinMCF final : public MinimumCostFlow {
 public:
    // flow
    OrlinMCF(NetworKit::Graph const& g, edgeid_map<long long> const& costs,
        NetworKit::node s, NetworKit::node t, long long f) : MinimumCostFlow(g, costs, s, t, f) {
        initFlow();
    }
    // flow
    OrlinMCF(NetworKit::Graph const& g, edgeid_map<long long> const& costs,
            node_map<long long> const& b) : MinimumCostFlow(g, costs, b) {
        initFlow();
    }
    // circulation
    OrlinMCF(NetworKit::Graph const& g, edgeid_map<std::pair<long long,long long>> const& bounds,
            edgeid_map<long long> const& costs) : MinimumCostFlow(g, bounds, costs) {
        initCirculation();
    }

 private:
    void makeCostsNonnegative();

    const double ALPHA = 0.9;
    long long delta = 0;

    void runImpl() override;
    void initFlow() override;
    void initCirculation() override;
    long long cp(NetworKit::node, NetworKit::node, NetworKit::edgeid);

    void augment(NetworKit::node, NetworKit::node, NetworKit::edgeid, bool, long long);
    void augment(NetworKit::node, NetworKit::node, 
        const std::vector<std::pair<NetworKit::edgeid, bool>>&, long long);
    
    NetworKit::Graph generateGraphForSP();
    bool isImbalanced();
    void initialize();

    std::stack<std::pair<NetworKit::node, NetworKit::node>> contractions;
    void contractNodes(NetworKit::node, NetworKit::node);
    void deltaScalingPhase();
    void contractionPhase();
    void uncontractAndComputeFlow();
    void scalingAugmentation(NetworKit::node, NetworKit::node);
    void fillDistancesUncapacitated(std::vector<double>&);
    NetworKit::Graph originalGraph;
};

} /* namespace Koala */

 