#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <algorithm>

namespace Koala {

class EdmondsKarpMCF final : public MinimumCostFlow {
 public: 
   using MinimumCostFlow::MinimumCostFlow;

 private:
   void initFlow() override {}
   void initCirculation() override {
      constructFlowFromCirculation();
   };

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