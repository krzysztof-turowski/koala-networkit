#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <flow/MinimumCostFlow.hpp>

namespace Koala {

class EdmondsKarpMinimumCostFlow final : public MinimumCostFlow {
  struct Edge {
      NetworKit::node from, to;
      int64_t cost, capacity, flow;
  };
  std::vector<Edge> edges;
  std::vector<std::vector<NetworKit::index>> neighbors;
  std::vector<int64_t> b, excess;
  std::vector<int64_t> potential;
  std::unordered_map<NetworKit::Edge, std::int64_t> computed_flow;
  NetworKit::count n;

  void run_impl() override;
  void initialize();
  void delta_scaling_phase(int64_t);
  void augmenting_phase(NetworKit::node, NetworKit::node, int64_t);
  void send(NetworKit::index, int64_t);
  std::vector<std::pair<int64_t, NetworKit::index>> dijkstra(
      NetworKit::node source, int64_t delta);

 public:
  EdmondsKarpMinimumCostFlow(MCFlowNetwork const& network) : MinimumCostFlow(network) {}
  int64_t getFlow(NetworKit::Edge const&) override;
  std::unordered_map<NetworKit::Edge, std::int64_t> getMinCostFlow() const override {
      return computed_flow;
  }
};

} /* namespace Koala */
