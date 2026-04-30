#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>

#include "io/DimacsGraphReader.hpp"
#include "graph/GraphTools.hpp"

#include "flow/electrical_flow/ElectricalFlow.hpp"

#include "../helpers.hpp"

class GenTest : public testing::Test {};

TEST(GenTest, testSuccess) {
  auto [G, s, t] = Koala::DimacsGraphReader().read_all("input/example.flow");
  G = Koala::GraphTools::convertDirectedGraphToUndirected(G, true);
  G = Koala::GraphTools::convertUndirectedGraphToDirected(G, true);
  int F = 15;

  Koala::ElectricalFlow ef(G, s, t);
  ef.run();

  auto graph = ef.getGraph();
  auto flow = ef.getFlow();

  int N = graph.numberOfNodes();
  for (int u = 0; u < N; ++u) {
    double demand = 0;
    for (int v = 0; v < N; ++v) {
      demand += flow[u][v];
      EXPECT_LE(abs(flow[u][v]), graph.weight(u, v));
      EXPECT_EQ(flow[u][v], -flow[v][u]);
      std::cout << flow[u][v] << "/" << graph.weight(u, v) << "\t\t";
    }
    std::cout << '\n';

    if (u == s) {
      EXPECT_EQ(demand, -F);
    } else if (u == t) {
      EXPECT_EQ(demand, F);
    } else {
      EXPECT_EQ(demand, 0);
    }
  }

  EXPECT_EQ(ef.getFlowSize(), F);
}
