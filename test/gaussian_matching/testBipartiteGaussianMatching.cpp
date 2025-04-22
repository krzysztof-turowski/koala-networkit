#include <gtest/gtest.h>
#include <vector>

#include <matching/BipartiteGaussianMatching.hpp>

#include "../helpers.hpp"

using namespace std;

// Example graph definition
const int N = 10;
const list<pair<int, int>> adj = {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 4},
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 8},
    {8, 9},
};

class BpTest : public testing::Test {};

TEST(BpTest, testSuccess) {
    auto G = build_graph(N, adj, false);
    auto bp = Koala::BipartiteGaussianMatching(G);

    bp.run();

    auto M = bp.getMatching();
    EXPECT_EQ(M.size(), N / 2);

    vector<int> counts(N, 0);
    for (auto m : M) {
        counts[m.first]++;
        counts[m.second]++;
    }

    for (auto c : counts) {
        EXPECT_EQ(c, 1);
    }
}

TEST(BpTest, testFail) {
    auto G = build_graph(N, adj, false);
    G.removeEdge(0, 1);

    auto bp = Koala::BipartiteGaussianMatching(G);

    bp.run();
    auto M = bp.getMatching();
    EXPECT_EQ(M.size(), 0);
}
