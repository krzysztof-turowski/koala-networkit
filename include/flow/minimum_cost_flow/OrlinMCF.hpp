#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Attributes.hpp>

namespace Koala {

class OrlinMCF : public MinimumCostFlow {
 public:
    
    OrlinMCF(const NetworKit::Graph&, const edgeid_map<int>&, const node_map<int>&);

 private:
    void makeCostsNonnegative();

    const double ALPHA = 0.9;
    int delta = 0;

    void runImpl() override;
    int cp(NetworKit::node, NetworKit::node, NetworKit::edgeid);

    void augment(NetworKit::node, NetworKit::node, NetworKit::edgeid, bool, int);
    void augment(NetworKit::node, NetworKit::node, 
        const std::vector<std::pair<NetworKit::edgeid, bool>>&, int);
    
    NetworKit::Graph generateGraphForSP();
    bool isImbalanced();
    void initialize();

    std::stack<std::pair<NetworKit::node, NetworKit::node>> contractions;
    void contractNodes(NetworKit::node, NetworKit::node);
    void deltaScalingPhase();
    void contractionPhase();
    void uncontractAndComputeFlow();
    void scalingAugmentation(NetworKit::node, NetworKit::node);
    edgeid_map<int> potential;
    NetworKit::Graph originalGraph;
    int delta = 0;
};

} /* namespace Koala */

 