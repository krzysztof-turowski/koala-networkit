#include <gtest/gtest.h>

#include <list>

#include "pathwidth/CographPathwidth.hpp"

#include "helpers.hpp"

struct MaxCliqueParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template<typename Algorithm>
class SimpleGraphsPathwidth : public testing::Test {
 protected:
    void verify(MaxCliqueParameters &parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E);
        NetworKit::count pathwidth;
        if constexpr (std::is_same_v<Algorithm, Koala::CographPathwidth>) {
            auto recognition = Koala::HabibPaulCographRecognition(G);
            recognition.run();
            if (recognition.isCograph()) {
                auto algorithm = Algorithm(G, recognition.cotree);
                algorithm.run();
                pathwidth = algorithm.getPathwidthSize();
            }
        } else {
            auto algorithm = Algorithm(G);
            algorithm.run();
            pathwidth = algorithm.getPathwidthSize();
        }
        EXPECT_EQ(pathwidth, parameters.expectedSetSize);
    }
};

TYPED_TEST_CASE_P(SimpleGraphsPathwidth);

TYPED_TEST_P(SimpleGraphsPathwidth, Cograph) {
    MaxCliqueParameters parameters =
        {6, {{0, 2}, {0, 3}, {0, 4}, {0, 1}, {5, 1}, {5, 2}, {5, 3}, {5, 4}, {1, 0}, {2, 0},
             {3, 0}, {4, 0}, {1, 5}, {2, 5}, {3, 5}, {4, 5}}, 8};
    this->verify(parameters);
}

REGISTER_TYPED_TEST_CASE_P(SimpleGraphsPathwidth, Cograph);

using Algorithms = testing::Types<
    Koala::CographPathwidth
>;

INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphsPathwidth, Algorithms);
