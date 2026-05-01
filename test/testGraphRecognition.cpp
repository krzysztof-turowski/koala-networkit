#include <gtest/gtest.h>

#include <list>

#include <recognition/CographRecognition.hpp>
#include <recognition/CographRecognitionOther.hpp>
#include <recognition/PerfectGraphRecognition.hpp>
#include <recognition/ChordalGraphRecognition.hpp>

#include "helpers.hpp"

struct GraphRecognitionParameters {
    int N;
    std::list<std::pair<int, int>> E;
    bool is_recognized;
};

auto perfect_graphs = testing::Values(
    GraphRecognitionParameters{4, {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {3, 1}, {4, 3}}, true},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {2, 4}, {4, 5}, {5, 6}, {2, 6}, {2, 3}, {3, 4},
                                   {0, 2}, {1, 4}}, true}
);

class PerfectGraphRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(PerfectGraphRecognitionTest, test_perfect_graphs) {
    auto parameters = GetParam();
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Koala::PerfectGraphRecognition(G);
    algorithm.run();

    auto is_perfect = algorithm.isPerfect();
    EXPECT_EQ(is_perfect, parameters.is_recognized);
}

INSTANTIATE_TEST_SUITE_P(
    test_perfect_graphs_example, PerfectGraphRecognitionTest, perfect_graphs);

auto cographs = testing::Values(
    GraphRecognitionParameters{4, {{0, 2}, {0, 3}, {1, 3}, {2, 3}}, true},
    GraphRecognitionParameters{5, {{0, 3}, {0, 4}, {3, 4}, {1, 3}, {2, 4}}, false},
    GraphRecognitionParameters{3, {{0, 1}, {0, 2}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 4}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {3, 4}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {3, 4}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {1, 4}, {2, 4}}, false},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                   {4, 5}, {5, 6}, {2, 6}}, false},
    GraphRecognitionParameters{7, {{0, 1}, {1, 2}, {0, 2}, {1, 4}, {2, 4}, {3, 2}, {3, 4},
                                   {4, 5}, {5, 6}}, false}
);

template <typename Algorithm>
void CographRecognitionTest(GraphRecognitionParameters parameters) {
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Algorithm(G);
    algorithm.run();

    auto is_cograph = algorithm.isCograph();
    EXPECT_EQ(is_cograph, parameters.is_recognized);
}

class CographBCHPRecognitionTest : public testing::TestWithParam<GraphRecognitionParameters> { };
class CographCSPRecognitionTest : public testing::TestWithParam<GraphRecognitionParameters> { };
class CographDahlhausRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };
class CographHabibPaulRecognitionTest
    : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(CographBCHPRecognitionTest, test_cographs) {
    CographRecognitionTest<Koala::BretscherCorneilHabibPaulCographRecognition>(GetParam());
}

TEST_P(CographCSPRecognitionTest, test_cographs) {
    CographRecognitionTest<Koala::CorneilStewartPerlCographRecognition>(GetParam());
}

TEST_P(CographDahlhausRecognitionTest, test_cographs) {
    CographRecognitionTest<Koala::DahlhausCographRecognition>(GetParam());
}

TEST_P(CographHabibPaulRecognitionTest, test_cographs) {
    CographRecognitionTest<Koala::HabibPaulCographRecognition>(GetParam());
}

INSTANTIATE_TEST_SUITE_P(test_cographs_bchp_example, CographBCHPRecognitionTest, cographs);
INSTANTIATE_TEST_SUITE_P(test_cographs_csp_example, CographCSPRecognitionTest, cographs);
INSTANTIATE_TEST_SUITE_P(test_cographs_dahlhaus_example, CographDahlhausRecognitionTest, cographs);
INSTANTIATE_TEST_SUITE_P(
    test_cographs_habib_paul_example, CographHabibPaulRecognitionTest, cographs);

auto chordal_graphs = testing::Values(
    GraphRecognitionParameters{4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}, true},
    GraphRecognitionParameters{4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, false},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, false},
    GraphRecognitionParameters{4, {{0, 1}, {0, 2}, {0, 3}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, 
                                   {1, 2}, {1, 3}, {1, 4}, 
                                   {2, 3}, {2, 4}, 
                                   {3, 4}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 0}, {3, 4}}, true},
    GraphRecognitionParameters{6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, 
                                   {1, 4}, {1, 5}, {2, 4}}, true},
    GraphRecognitionParameters{5, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}}, false}
);

template <typename Algorithm>
void ChordalRecognitionTest(GraphRecognitionParameters parameters) {
    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);
    auto algorithm = Algorithm(G);
    algorithm.run();

    auto is_chordal = algorithm.isChordal();
    EXPECT_EQ(is_chordal, parameters.is_recognized);
}

class ChordalGraphMCSRecognitionTest : public testing::TestWithParam<GraphRecognitionParameters> { };
class ChordalGraphLexBFSRecognitionTest : public testing::TestWithParam<GraphRecognitionParameters> { };

TEST_P(ChordalGraphMCSRecognitionTest, test_chordal_graphs) {
    ChordalRecognitionTest<Koala::MaximumCardinalitySearchChordalGraphRecognition>(GetParam());
}

TEST_P(ChordalGraphLexBFSRecognitionTest, test_chordal_graphs) {
    ChordalRecognitionTest<Koala::LexBFSChordalGraphRecognition>(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
    test_chordal_graphs_mcs, 
    ChordalGraphMCSRecognitionTest, 
    chordal_graphs
);

INSTANTIATE_TEST_SUITE_P(
    test_chordal_graphs_lexbfs, 
    ChordalGraphLexBFSRecognitionTest, 
    chordal_graphs
);
