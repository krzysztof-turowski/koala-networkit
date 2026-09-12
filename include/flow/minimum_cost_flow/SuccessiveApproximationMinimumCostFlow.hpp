#pragma once

#include <cstdint>
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <flow/MinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/MCFlowNetwork.hpp>

namespace Koala {

class SuccessiveApproximationMinimumCostFlow final : public MinimumCostFlow {
  struct Edge {
    NetworKit::node from, to;
    int64_t cost, capacity, flow;
  };
  std::vector<Edge> edges;
  std::vector<std::vector<NetworKit::index>> neighbors;
  std::unordered_map<NetworKit::Edge, int64_t> computed_flow;
  void run_impl() override;
  bool is_imbalanced();
  void initialize();
  void push(NetworKit::index edge_index);
  void relabel(NetworKit::node const& u);
  void refine();
  void wave();
  bool discharge(NetworKit::node const& u);

  double reduced_cost(NetworKit::index edge_index);
  int64_t residual_capacity(NetworKit::index edge_index);
  void force_flow(NetworKit::index edge_index, int64_t amount);
  std::vector<double> potential;
  std::vector<int64_t> excess;
  NetworKit::count nodes_number{0};

  double epsilon{0.};

  class DischargeList {
   public:
    virtual NetworKit::node getNext() { return 0; }
    virtual void moveToStart() {}
    virtual ~DischargeList() = default;
  };

  class ToposortList : public DischargeList {
   public:
    explicit ToposortList(SuccessiveApproximationMinimumCostFlow&);

    NetworKit::node getNext() override;
    void moveToStart() override;

   private:
    SuccessiveApproximationMinimumCostFlow &algorithm;
    std::list<NetworKit::node> nodes;

    std::list<NetworKit::node>::iterator current;
    std::list<NetworKit::node>::iterator next;
  };

 public:
  explicit SuccessiveApproximationMinimumCostFlow(MCFlowNetwork const& network)
      : MinimumCostFlow(network) {}
  int64_t getFlow(NetworKit::Edge const& edge) override;
  std::unordered_map<NetworKit::Edge, int64_t> getMinCostFlow() const override;
};

}  /* namespace Koala */
