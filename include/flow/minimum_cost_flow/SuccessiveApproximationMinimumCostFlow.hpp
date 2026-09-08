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
  void push(NetworKit::index eid);
  void relabel(NetworKit::node const&);
  void refine();
  void wave();
  bool discharge(NetworKit::node const&);

  double cp(NetworKit::index eid);
  int64_t uf(NetworKit::index eid);
  void force_flow(NetworKit::index eid, int64_t f);
  std::vector<double> potential;
  std::vector<int64_t> excess;
  NetworKit::count nodes_number{0};

  double epsi{0.};

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
    SuccessiveApproximationMinimumCostFlow &approx;
    std::list<NetworKit::node> nodes;
    std::vector<bool> vis;
    void dfs(NetworKit::node);

    std::list<NetworKit::node>::iterator it1;
    std::list<NetworKit::node>::iterator it2;
  };

 public:
  explicit SuccessiveApproximationMinimumCostFlow(const MCFlowNetwork& network) : MinimumCostFlow(network) {}
  int64_t getFlow(const NetworKit::Edge& edge) override;
  std::unordered_map<NetworKit::Edge, std::int64_t> getMinCostFlow() const override {
      return computed_flow;
  }
};

}  /* namespace Koala */
