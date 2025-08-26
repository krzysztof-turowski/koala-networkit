#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <algorithm>

namespace Koala {

class EdmondsKarpMCF : public MinimumCostFlow {
 private:
    void runImpl() override;
    void initialize();
    void deltaScalingPhase();
    void send(NetworKit::node, NetworKit::node, NetworKit::edgeid, int);
    NetworKit::Graph getDeltaResidual();
    int cp(NetworKit::node, NetworKit::node, NetworKit::edgeid);

    edgeid_map<int> potential; 
    int delta = 0;


};

} /* namespace Koala */