#include <gtest/gtest.h>
#include <list>

#include "planar_sssp/PlanarSSSP.hpp"
#include "planar_sssp/FredericksonPlanarSSSP.hpp"
#include "helpers.hpp"

struct SSSPParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
    NetworKit::node s;
    NetworKit::node t;
    int expectedDistance;
};

class PlanarSSSPTest : public testing::TestWithParam<SSSPParameters> {
};

TEST_P(PlanarSSSPTest, TestPlanarSSSPSolution) {
    SSSPParameters const& parameters = GetParam();

    std::cout << "Test name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);

    Koala::FredericksonPlanarSSSP algorithm(G, parameters.s, parameters.t);

    algorithm.run();

    EXPECT_EQ(algorithm.getSSSDistanceToTarget(), parameters.expectedDistance);
}

INSTANTIATE_TEST_SUITE_P(
    Default, PlanarSSSPTest, testing::Values(SSSPParameters{ "star", 5, {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 6}}, 0, 2, 5 }
        , SSSPParameters{ "double-star", 8, {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 6}, {4,7,1}, {4,6,1}, { 4,5,1 }}, 0, 2, 5 }


        // SSSPParameters{ "", 4, {{0, 1, 1}, {0, 3, 2}, {1, 2, 3}, {2, 3, 4}}, 0, 2, 5 },

        // SSSPParameters{ "", 4, {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {1, 2, 1}, {1, 3, 1}, {2, 3, 1}}, 0, 2, 5 },

        // SSSPParameters{ "", 5, {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}, {0, 4, 1}}, 0, 2, 5 },

        // SSSPParameters{ "", 4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}}, 0, 2, 5 },

        // SSSPParameters{ "", 5, {{0, 1, 10}, {0, 4, 3}, {1, 2, 2}, {2, 3, 4}, {3, 4, 5}}, 0, 2, 5 }
    ));
    