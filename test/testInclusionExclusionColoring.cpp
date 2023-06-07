#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <coloring/InclusionExclusionVertexColoring.hpp>

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    std::vector<int> subset;
    int numberOfSets;
};

struct InclusionExclusionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int k;
};

class IndependentSetCheckerTest : public ::testing::TestWithParam<IndependentSetParameters> {};
class InclusionExclusionVertexColoringTest
: public ::testing::TestWithParam<InclusionExclusionParameters> {};

NetworKit::Graph build_graph(const int& N, const std::list<std::pair<int, int>>& E) {
    NetworKit::Graph G(N, false, false);
    for (const auto& [u, v] : E) {
        G.addEdge(u, v);
    }
    return G;
}


TEST_P(IndependentSetCheckerTest, testNumberOfIndependentSetsNotIntersectingWith) {
    auto parameters = GetParam();
    auto G = build_graph(parameters.N, parameters.E);
    Koala::IndependentSetChecker independentSetChecker(G);
    std::vector<NetworKit::node> subset(parameters.subset.begin(), parameters.subset.end());
    auto result = independentSetChecker.numberOfIndependentSetsNotIntersectingWith(subset);
    EXPECT_EQ(result, parameters.numberOfSets);
}

INSTANTIATE_TEST_SUITE_P(test,
IndependentSetCheckerTest,
testing::Values(IndependentSetParameters{ 7,
{ { 0, 1 }, { 0, 2 }, { 1, 3 }, { 1, 4 }, { 2, 3 }, { 2, 5 }, { 2, 6 }, { 3, 4 }, { 4, 5 },
{ 5, 6 } },
{ 2, 3 }, 13 }));

TEST_P(InclusionExclusionVertexColoringTest, testGetColoring) {
    auto parameters = GetParam();
    auto G = build_graph(parameters.N, parameters.E);
    Koala::InclusionExclusionVertexColoring coloring(G);
    coloring.run();
    auto result = coloring.getChromaticNumber();
    EXPECT_EQ(result, parameters.k);
}

INSTANTIATE_TEST_SUITE_P(test,
InclusionExclusionVertexColoringTest,
testing::Values(InclusionExclusionParameters{ 4, { { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 } }, 2 },
InclusionExclusionParameters{
6, { { 0, 1 }, { 0, 2 }, { 1, 2 }, { 0, 3 }, { 3, 4 }, { 1, 5 }, { 4, 5 } }, 3 },
InclusionExclusionParameters{ 10,
{ { 0, 2 }, { 0, 3 }, { 0, 4 }, { 0, 5 }, { 1, 2 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 2, 6 },
{ 2, 8 }, { 3, 4 }, { 3, 6 }, { 4, 6 }, { 5, 6 }, { 5, 7 }, { 5, 9 }, { 7, 8 }, { 7, 9 },
{ 8, 9 } },
3 },
InclusionExclusionParameters{ 8,
{ { 0, 1 }, { 0, 2 }, { 0, 4 }, { 0, 6 }, { 1, 2 }, { 1, 3 }, { 1, 7 }, { 2, 3 }, { 2, 4 },
{ 3, 5 }, { 3, 7 }, { 4, 5 }, { 4, 6 }, { 5, 6 }, { 5, 7 }, { 6, 7 } },
4 },
InclusionExclusionParameters{ 8,
{ { 0, 1 }, { 0, 2 }, { 0, 4 }, { 0, 6 }, { 1, 3 }, { 1, 5 }, { 1, 7 }, { 2, 3 }, { 2, 4 },
{ 2, 5 }, { 2, 6 }, { 3, 4 }, { 3, 5 }, { 3, 7 }, { 4, 6 }, { 4, 7 }, { 5, 6 }, { 5, 7 },
{ 6, 7 } },
4 },
InclusionExclusionParameters{ 10,
{ { 0, 1 }, { 0, 2 }, { 0, 3 }, { 0, 5 }, { 0, 6 }, { 0, 7 }, { 1, 2 }, { 1, 3 }, { 1, 4 },
{ 1, 6 }, { 1, 7 }, { 2, 3 }, { 2, 4 }, { 2, 5 }, { 2, 7 }, { 3, 4 }, { 3, 5 }, { 3, 6 }, { 4, 5 },
{ 4, 7 }, { 4, 8 }, { 4, 9 }, { 5, 6 }, { 5, 8 }, { 5, 9 }, { 6, 7 }, { 6, 8 }, { 6, 9 }, { 7, 8 },
{ 7, 9 }, { 8, 9 } },
5 },
InclusionExclusionParameters{ 4, { { 0, 2 }, { 1, 3 }, { 2, 3 } }, 2 },
InclusionExclusionParameters{ 9,
{ { 0, 4 }, { 0, 5 }, { 0, 6 }, { 0, 8 }, { 1, 5 }, { 1, 6 }, { 1, 7 }, { 1, 8 }, { 2, 6 },
{ 3, 7 }, { 3, 8 }, { 4, 5 }, { 4, 7 }, { 4, 8 }, { 5, 7 }, { 5, 8 }, { 6, 7 }, { 7, 8 } },
4 }));