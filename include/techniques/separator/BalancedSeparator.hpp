#pragma once

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

class BalancedSeparator : public NetworKit::Algorithm {
 public:
  explicit BalancedSeparator(const NetworKit::Graph& graph);

  struct Partition {
    std::vector<NetworKit::node> separator;
    std::vector<NetworKit::node> A;
    std::vector<NetworKit::node> B;
  };

  void run() override = 0;

  const Partition& getPartition() const;

  const std::vector<NetworKit::node>& getSeparator() const;
  const std::vector<NetworKit::node>& getPartitionA() const;
  const std::vector<NetworKit::node>& getPartitionB() const;

 protected:
  const NetworKit::Graph& graph;
  Partition partition;
};
