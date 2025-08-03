#include <gtest/gtest.h>
#include <vector>
#include <random>

#include <matching/gaussian_matching/GeneralGaussianMatching.hpp>
#include <io/DimacsGraphReader.hpp>

#include "../helpers.hpp"

using namespace std;

class GenTest : public testing::Test { };

TEST(GenTest, testSuccess) {
    srand(time(NULL));
    auto G = Koala::DimacsGraphReader().read("input/perfect.in");
    int n = G.numberOfNodes();

    Koala::GeneralGaussianMatching gen(G);
    gen.run();

    auto M = gen.getMatching();
    EXPECT_EQ(M.size(), n / 2);

    vector<int> counts(n, 0);
    for (auto [u,v] : M) {
        counts[u]++;
        counts[v]++;
        EXPECT_TRUE(G.hasEdge(u,v));
    }

    for (auto c : counts) {
        EXPECT_EQ(c, 1);
    }
}
