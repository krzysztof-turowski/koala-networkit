#pragma once

#include <climits>
#include <cstdint>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include <flow/GoldbergTarjanPushRelabelMaximumFlow.hpp>
#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Attributes.hpp>

namespace Koala {

class OrlinMinimumCostFlow final : public MinimumCostFlow {
  struct Edge {
    NetworKit::node from, to;
    int64_t cost, capacity, flow;
  };

  NetworKit::Graph original_graph;
  std::vector<Edge> edges;
  std::vector<Edge> original_edges;
  std::vector<std::vector<NetworKit::index>> neighbors;
  std::vector<int64_t> potential;
  // potential that is computed for the original network
  std::vector<int64_t> potential_computed;
  std::vector<int64_t> excess;
  std::vector<std::pair<int64_t, NetworKit::index>> dist;
  std::vector<bool> visited;
  NetworKit::count nodes_number;
  NetworKit::count max_nodeid;
  // Storing with costs for potential resotration
  std::stack<std::tuple<NetworKit::node, NetworKit::node, int64_t>> contracted_nodes;
  const double ALPHA = 0.7;

  std::optional<Koala::GoldbergTarjanPushRelabelMaximumFlow> maxflow;
  void run_impl() override;
  bool is_imbalanced();
  void initialize();
  int64_t find_optimal_delta(int64_t delta);
  void push_no_excess(NetworKit::index edge_idx, int64_t amount);
  void contract_nodes(NetworKit::node u, NetworKit::node v);
  void apply_potential();
  void contraction_phase(int64_t delta);
  void augmenting_phase(NetworKit::node s, NetworKit::node t, int64_t delta);
  void uncontract_nodes_potential();
  void dijkstra(NetworKit::node source, int64_t delta);
  void make_reduced_costs_nonnegative();
  std::stack<std::pair<NetworKit::node, NetworKit::node>> contractions;
  std::pair<NetworKit::node, NetworKit::node> uncapacitated_nodes_bounds;
  bool is_added_uncapacitated(NetworKit::node v) const;
  void compute_final_flows();
  std::unordered_map<NetworKit::Edge, std::int64_t> computed_flow;

 public:
  explicit OrlinMinimumCostFlow(Koala::MCFlowNetwork& network) : MinimumCostFlow(network) {}
  int64_t getFlow(const NetworKit::Edge& edge) override;
  std::unordered_map<NetworKit::Edge, std::int64_t> getMinCostFlow() const override {
      return computed_flow;
  }
};

}  /* namespace Koala */
