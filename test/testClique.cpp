#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <clique/ChordalClique.hpp>
#include <clique/CographClique.hpp>
#include <recognition/ChordalGraphRecognition.hpp>
#include <recognition/CographRecognition.hpp>

#include "test/helpers.hpp"

struct MaxCliqueParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template<typename Algorithm>
class SimpleGraphsClique : public testing::Test {
 protected:
    void verify(MaxCliqueParameters &parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E);
        std::vector<NetworKit::node> max_clique;
        if constexpr (std::is_same_v<Algorithm, Koala::CographMaxClique>) {
            auto recognition = Koala::HabibPaulCographRecognition(G);
            recognition.run();
            if (recognition.isCograph()) {
                auto algorithm = Algorithm(G, recognition.cotree);
                algorithm.run();
                max_clique = algorithm.getMaxClique();
            }
        } else {
            auto algorithm = Algorithm(G);
            algorithm.run();
            max_clique = algorithm.getMaxClique();
        }
        for (auto x : max_clique) {
            for (auto y : max_clique) {
                EXPECT_TRUE((x == y) || (G.hasEdge(x, y)));
            }
        }
        EXPECT_EQ(max_clique.size(), parameters.expectedSetSize);
    }
};

TYPED_TEST_CASE_P(SimpleGraphsClique);

TYPED_TEST_P(SimpleGraphsClique, Cograph) {
    MaxCliqueParameters parameters = {
        6, {{0, 2}, {0, 3}, {0, 4}, {0, 1}, {5, 1}, {5, 2}, {5, 3}, {5, 4}, {1, 0}, {2, 0},
            {3, 0}, {4, 0}, {1, 5}, {2, 5}, {3, 5}, {4, 5}}, 2};
    this->verify(parameters);
}

REGISTER_TYPED_TEST_CASE_P(SimpleGraphsClique, Cograph);

using Algorithms = testing::Types<Koala::CographMaxClique>;

INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphsClique, Algorithms);

class ChordalMaxCliqueTest
    : public testing::TestWithParam<MaxCliqueParameters> { };

TEST_P(ChordalMaxCliqueTest, test) {
    MaxCliqueParameters const &parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto recognition = Koala::MaximumCardinalitySearchChordalGraphRecognition(G);
    recognition.run();
    EXPECT_TRUE(recognition.isChordal());
    auto algorithm = Koala::ChordalMaxClique(G, recognition.getPEO());
    algorithm.run();
    auto clique = algorithm.getMaxClique();
    for (auto x : clique)
        for (auto y : clique)
            EXPECT_TRUE(x == y || G.hasEdge(x, y));
    EXPECT_EQ(static_cast<int>(clique.size()), parameters.expectedSetSize);
}

INSTANTIATE_TEST_SUITE_P(
    test_example, ChordalMaxCliqueTest, testing::Values(
    MaxCliqueParameters{0, {}, 0},
    MaxCliqueParameters{1, {}, 1},
    MaxCliqueParameters{4, {{0, 1}, {0, 2}, {0, 3}}, 2},
    MaxCliqueParameters{4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}, 3},
    MaxCliqueParameters{
        5, {{0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 5},
    MaxCliqueParameters{5, {{0, 1}, {1, 2}, {2, 0}, {3, 4}}, 3},
    MaxCliqueParameters{
        6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
            {1, 4}, {1, 5}, {2, 4}}, 3}
));
