#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <algorithm>

namespace Koala {

class EdmondsKarpMCF final : public MinimumCostFlow {
 public: 
  // flow
  EdmondsKarpMCF(NetworKit::Graph const& g, edge_map<long long> const& costs,
      NetworKit::node s, NetworKit::node t, long long f) : MinimumCostFlow(g, costs, s, t, f) {
      initFlow();
  }
  // flow
  EdmondsKarpMCF(NetworKit::Graph const& g, edge_map<long long> const& costs,
          node_map<long long> const& b) : MinimumCostFlow(g, costs, b) {
      initFlow();
  }
  // circulation
  EdmondsKarpMCF(NetworKit::Graph const& g, edge_map<std::pair<long long,long long>> const& bounds,
          edge_map<long long> const& costs) : MinimumCostFlow(g, bounds, costs) {
      initCirculation();
  }

 private:
   void initFlow() override {}
   void initCirculation() override {
      constructFlowFromCirculation();
   };

   void runImpl() override;
   void initialize();
   void deltaScalingPhase();
   void send(NetworKit::node, NetworKit::node, long long);
   NetworKit::Graph getDeltaResidual();
   long long cp(NetworKit::node, NetworKit::node);

   long long delta = 0;
};

} /* namespace Koala */