#include <gtest/gtest.h>

#include <list>

#include "max_clique/CographMaxClique.hpp"
#include "helpers.hpp"

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
        std::set<NetworKit::node> max_clique;
        if constexpr (std::is_same_v<Algorithm, Koala::CographMaxClique>) {
            auto recognition = Koala::CographRecognition(G);
            recognition.run();
            if (recognition.isCograph()) {
                auto algorithm = Algorithm(G, recognition.cotree);
                algorithm.run();
                max_clique = algorithm.getMaxCliqueSet();
            }
        } else {
            auto algorithm = Algorithm(G);
            algorithm.run();
            max_clique = algorithm.getMaxCliqueSet();
        }
        for (auto x : max_clique) {
            for (auto y : max_clique) {
                EXPECT_TRUE((x == y) || (G.hasEdge(x, y)));
            }
        }
        EXPECT_EQ(max_clique.size(), parameters.expectedSetSize);
    }
};

TYPED_TEST_CASE_P
(SimpleGraphsClique);

TYPED_TEST_P(SimpleGraphsClique, Cograph) {
    MaxCliqueParameters parameters =
            {6, {
                    {0, 2}, {0, 3}, {0, 4}, {0, 1}, {5, 1}, {5, 2}, {5, 3}, {5, 4}, {1, 0}, {2, 0}, {3, 0}, {4, 0},
                    {1, 5}, {2, 5}, {3, 5}, {4, 5}
            },
             2};
    this->verify(parameters);
}

REGISTER_TYPED_TEST_CASE_P
(
        SimpleGraphsClique, Cograph);

using Algorithms = testing::Types<
        Koala::CographMaxClique
>;

INSTANTIATE_TYPED_TEST_CASE_P
(IndependentSet, SimpleGraphsClique, Algorithms);
