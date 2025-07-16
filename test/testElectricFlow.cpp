#include <gtest/gtest.h>

#include <electric_flow/ElectricFlow.hpp>
#include <io/DimacsGraphReader.hpp>

#include "helpers.hpp"

using namespace std;

class GenTest : public testing::Test {};

TEST(GenTest, testSuccess) {
  auto G = Koala::DimacsGraphReader().read("input/flow.dat");

  Koala::ElectricFlow ef(G, 5, 4);
  ef.run();

  int N = ef.graph.numberOfNodes();
  for (int u = 0; u < N; ++u) {
    for (int v = 0; v < N; ++v) {
      cout << round(ef.primal.flow[u][v]*100.0)/100.0 << "/" << ef.graph.weight(u, v) << "\t\t";
    }
    cout << '\n';
  }

  EXPECT_EQ(ef.getMaxFlow(), 100);
}
