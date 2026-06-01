#pragma once

#include <utility>
#include <vector>

#include <structures/DynamicTree.hpp>

#include <networkit/graph/Graph.hpp>

namespace Koala {

template <typename Value>
class NaiveDynamicTree final : public DynamicTree<Value> {
 public:
  NaiveDynamicTree(NetworKit::count n, std::vector<std::vector<Value>> &weights);

  void link(NetworKit::node u, NetworKit::node v, Value value) override;
  void cut(NetworKit::node u, NetworKit::node v) override;
  NetworKit::node findRoot(NetworKit::node v) override;
  void pathAdd(NetworKit::node u, NetworKit::node v, Value value);
  NetworKit::Edge pathMin(NetworKit::node u, NetworKit::node v);
  Value pathSum(NetworKit::node u, NetworKit::node v);

  Value getWeight(NetworKit::node u, NetworKit::node v) const;
  void addWeight(NetworKit::node u, NetworKit::node v, Value value);

 private:
  NetworKit::Graph graph;
  std::vector<std::vector<Value>> &weights;
};

};  // namespace Koala
