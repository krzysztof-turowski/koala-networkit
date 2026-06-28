#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <clique_cover/ChordalCliqueCover.hpp>
#include <recognition/ChordalGraphRecognition.hpp>
#include <test/helpers.hpp>

struct CliqueCoverParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template<typename Algorithm>
class SimpleGraphsCliqueCover : public testing::Test {
 protected:
    void verify(CliqueCoverParameters &parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
        std::vector<std::vector<NetworKit::node>> cover;
        if constexpr (std::is_same_v<Algorithm, Koala::ChordalMinCliqueCover>) {
            auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
            recognition.run();
            if (recognition.isChordal()) {
                auto algorithm = Algorithm(G, recognition.getPEO());
                algorithm.run();
                cover = algorithm.getCliqueCover();
            }
        }
        std::vector<bool> covered(G.upperNodeIdBound(), false);
        for (const auto &part : cover) {
            for (auto x : part) {
                EXPECT_FALSE(covered[x]);
                covered[x] = true;
                for (auto y : part) EXPECT_TRUE(x == y || G.hasEdge(x, y));
            }
        }
        G.forNodes([&](NetworKit::node v) { EXPECT_TRUE(covered[v]); });
        EXPECT_EQ(static_cast<int>(cover.size()), parameters.expectedSetSize);
    }
};

TYPED_TEST_CASE_P(SimpleGraphsCliqueCover);

TYPED_TEST_P(SimpleGraphsCliqueCover, StarGraph) {
    CliqueCoverParameters parameters = {4, {{0, 1}, {0, 2}, {0, 3}}, 3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphsCliqueCover, CompleteGraph) {
    CliqueCoverParameters parameters = {
        5, {{0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 1};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphsCliqueCover, TriangleWithIsolatedEdge) {
    CliqueCoverParameters parameters = {5, {{0, 1}, {1, 2}, {2, 0}, {3, 4}}, 2};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphsCliqueCover, ChordalGraph) {
    CliqueCoverParameters parameters = {
        6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
            {1, 4}, {1, 5}, {2, 4}}, 2};
    this->verify(parameters);
}

REGISTER_TYPED_TEST_CASE_P(
    SimpleGraphsCliqueCover, StarGraph, CompleteGraph, TriangleWithIsolatedEdge, ChordalGraph);

using Algorithms = testing::Types<Koala::ChordalMinCliqueCover>;

INSTANTIATE_TYPED_TEST_CASE_P(CliqueCover, SimpleGraphsCliqueCover, Algorithms);

class ChordalMinCliqueCoverTest
    : public testing::TestWithParam<CliqueCoverParameters> { };

TEST_P(ChordalMinCliqueCoverTest, test) {
    CliqueCoverParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
    recognition.run();
    EXPECT_TRUE(recognition.isChordal());
    auto algorithm = Koala::ChordalMinCliqueCover(G, recognition.getPEO());
    algorithm.run();
    auto cover = algorithm.getCliqueCover();
    std::vector<bool> covered(G.upperNodeIdBound(), false);
    for (const auto &part : cover) {
        for (auto x : part) {
            EXPECT_FALSE(covered[x]);
            covered[x] = true;
            for (auto y : part) EXPECT_TRUE(x == y || G.hasEdge(x, y));
        }
    }
    G.forNodes([&](NetworKit::node v) { EXPECT_TRUE(covered[v]); });
    EXPECT_EQ(static_cast<int>(cover.size()), parameters.expectedSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, ChordalMinCliqueCoverTest, testing::Values(
    CliqueCoverParameters{0, {}, 0},
    CliqueCoverParameters{1, {}, 1},
    CliqueCoverParameters{4, {{0, 1}, {0, 2}, {0, 3}}, 3},
    CliqueCoverParameters{4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}, 2},
    CliqueCoverParameters{
        5, {{0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 1},
    CliqueCoverParameters{5, {{0, 1}, {1, 2}, {2, 0}, {3, 4}}, 2},
    CliqueCoverParameters{
        6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
            {1, 4}, {1, 5}, {2, 4}}, 2}
));
