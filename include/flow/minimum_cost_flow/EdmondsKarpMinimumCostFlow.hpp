#pragma once

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
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
  std::vector<int64_t> excess;
  std::vector<int64_t> potential;
  NetworKit::count max_node_id;

  void run_impl() override;
  void initialize();
  void delta_scaling_phase(int64_t delta);
  void augmenting_phase(NetworKit::node s, NetworKit::node t, int64_t delta);
  void send(NetworKit::index edge_index, int64_t amount);
  std::vector<std::pair<int64_t, NetworKit::index>> dijkstra(
      NetworKit::node source, int64_t delta);

 public:
  explicit EdmondsKarpMinimumCostFlow(MCFlowNetwork network) : MinimumCostFlow(network) {}
};

} /* namespace Koala */
