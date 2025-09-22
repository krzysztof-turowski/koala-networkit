#include <gtest/gtest.h>

#include <flow/electrical_flow/ElectricalFlow.hpp>
#include <io/DimacsGraphReader.hpp>
#include <networkit/graph/Graph.hpp>

#include "../helpers.hpp"

using namespace std;
using namespace NetworKit;

class GenTest : public testing::Test {
};

TEST(GenTest, testSuccess) {
  auto G = Koala::DimacsGraphReader().read("input/flow.dat");

  int s = 5, t = 4, F = 100;

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
      cout << flow[u][v] << "/" << graph.weight(u, v) << "\t\t";
    }
    cout << '\n';

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
