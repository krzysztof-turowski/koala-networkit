#include <gtest/gtest.h>

#include <list>

#include <io/D6GraphReader.hpp>
#include <io/D6GraphWriter.hpp>
#include <io/DimacsBinaryGraphReader.hpp>
#include <io/DimacsBinaryGraphWriter.hpp>
#include <io/DimacsGraphReader.hpp>
#include <io/DimacsGraphWriter.hpp>
#include <io/G6GraphReader.hpp>
#include <io/G6GraphWriter.hpp>
#include <io/S6GraphReader.hpp>
#include <io/S6GraphWriter.hpp>

#include "helpers.hpp"

struct GraphIOParameters {
    std::string G;
    int N;
    std::list<std::pair<int, int>> E;
};

class GraphReaderFromDigraph6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphReaderFromGraph6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphReaderFromSparse6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphReaderFromDimacsTest
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphReaderFromDimacsBinaryTest
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphWriterToDigraph6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphWriterToGraph6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphWriterToSparse6Test
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphWriterToDimacsTest
    : public testing::TestWithParam<GraphIOParameters> { };

class GraphWriterToDimacsBinaryTest
    : public testing::TestWithParam<GraphIOParameters> { };

TEST_P(GraphReaderFromDigraph6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G = Koala::D6GraphReader().readline(parameters.G);
    EXPECT_EQ(G.numberOfNodes(), parameters.N);
    EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
    for (const auto &e : parameters.E) {
        EXPECT_TRUE(G.hasEdge(e.first, e.second));
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphReaderFromDigraph6Test, testing::Values(
        GraphIOParameters{
            "&DI?AO?", 5, {{0, 2}, {0, 4}, {3, 1}, {3, 4}}}
));

TEST_P(GraphReaderFromGraph6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G = Koala::G6GraphReader().readline(parameters.G);
    EXPECT_EQ(G.numberOfNodes(), parameters.N);
    EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
    for (const auto &e : parameters.E) {
        EXPECT_TRUE(G.hasEdge(e.first, e.second));
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphReaderFromGraph6Test, testing::Values(
        GraphIOParameters{
            "G?r@`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "G?qa`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 2}, {6, 1}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "GCQR@O", 8, {{3, 0}, {4, 1}, {5, 0}, {5, 3}, {6, 1}, {6, 2}, {7, 2}, {7, 4}}},
        GraphIOParameters{
            "Fw??G", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}}
));

TEST_P(GraphReaderFromSparse6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G = Koala::S6GraphReader().readline(parameters.G);
    EXPECT_EQ(G.numberOfNodes(), parameters.N);
    EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
    for (const auto &e : parameters.E) {
        EXPECT_TRUE(G.hasEdge(e.first, e.second));
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphReaderFromSparse6Test, testing::Values(
        GraphIOParameters{
            ":Fa@x^", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}},
        GraphIOParameters{
            ":Go@_YMb", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}}
));

TEST_P(GraphReaderFromDimacsTest, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G = Koala::DimacsGraphReader().read(parameters.G);
    EXPECT_EQ(G.numberOfNodes(), parameters.N);
    EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
    for (const auto &e : parameters.E) {
        EXPECT_TRUE(G.hasEdge(e.first, e.second));
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphReaderFromDimacsTest, testing::Values(
        GraphIOParameters{
            "input/example_1.col", 8,
            {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "input/example_2.col", 11,
            {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 2}, {4, 3}, {5, 0}, {5, 1}, {5, 3}, {5, 4},
             {6, 0}, {6, 1}, {6, 2}, {7, 1}, {7, 2}, {7, 3}, {7, 5}, {7, 6}, {8, 2}, {8, 3},
             {8, 5}, {8, 6}, {9, 0}, {9, 2}, {9, 3}, {9, 5}, {9, 6}, {9, 8}, {10, 1},
             {10, 2}, {10, 3}, {10, 5}, {10, 6}}}
));

TEST_P(GraphReaderFromDimacsBinaryTest, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G = Koala::DimacsBinaryGraphReader().read(parameters.G);
    EXPECT_EQ(G.numberOfNodes(), parameters.N);
    EXPECT_EQ(G.numberOfEdges(), parameters.E.size());
    for (const auto &e : parameters.E) {
        EXPECT_TRUE(G.hasEdge(e.first, e.second));
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphReaderFromDimacsBinaryTest, testing::Values(
        GraphIOParameters{
            "input/example_1.col.b", 11,
            {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 2}, {4, 3}, {5, 0}, {5, 1}, {5, 3}, {5, 4},
             {6, 0}, {6, 1}, {6, 2}, {7, 1}, {7, 2}, {7, 3}, {7, 5}, {7, 6}, {8, 2}, {8, 3},
             {8, 5}, {8, 6}, {9, 0}, {9, 2}, {9, 3}, {9, 5}, {9, 6}, {9, 8}, {10, 1},
             {10, 2}, {10, 3}, {10, 5}, {10, 6}}}
));

TEST_P(GraphWriterToDigraph6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G(parameters.N, false, true);
    for (const auto &[u, v] : parameters.E) {
        G.addEdge(u, v);
    }
    EXPECT_EQ(Koala::D6GraphWriter().writeline(G), parameters.G);
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphWriterToDigraph6Test, testing::Values(
        GraphIOParameters{
            "&DI?AO?", 5, {{0, 2}, {0, 4}, {3, 1}, {3, 4}}}
));

TEST_P(GraphWriterToGraph6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G(parameters.N, false, false);
    for (const auto &[u, v] : parameters.E) {
        G.addEdge(u, v);
    }
    EXPECT_EQ(Koala::G6GraphWriter().writeline(G), parameters.G);
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphWriterToGraph6Test, testing::Values(
        GraphIOParameters{
            "G?r@`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "G?qa`_", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 2}, {6, 1}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "GCQR@O", 8, {{3, 0}, {4, 1}, {5, 0}, {5, 3}, {6, 1}, {6, 2}, {7, 2}, {7, 4}}},
        GraphIOParameters{
            "Fw??G", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}}
));

TEST_P(GraphWriterToSparse6Test, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G(parameters.N, false, false);
    for (const auto &[u, v] : parameters.E) {
        G.addEdge(u, v);
    }
    EXPECT_EQ(Koala::S6GraphWriter().writeline(G), parameters.G);
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphWriterToSparse6Test, testing::Values(
        GraphIOParameters{
            ":Fa@x^", 7, {{1, 0}, {2, 0}, {2, 1}, {6, 5}}},
        GraphIOParameters{
            ":Go@_YMb", 8, {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}}
));

TEST_P(GraphWriterToDimacsTest, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G(parameters.N, false, false);
    for (const auto &[u, v] : parameters.E) {
        G.addEdge(u, v);
    }
    std::string path = generate_filename("input");
    Koala::DimacsGraphWriter().write(G, path);
    EXPECT_TRUE(compare_files(path, parameters.G));
    remove(path.data());
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphWriterToDimacsTest, testing::Values(
        GraphIOParameters{
            "input/example_1.col", 8,
            {{4, 0}, {4, 1}, {5, 0}, {5, 1}, {6, 2}, {6, 3}, {7, 2}, {7, 3}}},
        GraphIOParameters{
            "input/example_2.col", 11,
            {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 2}, {4, 3}, {5, 0}, {5, 1}, {5, 3}, {5, 4},
             {6, 0}, {6, 1}, {6, 2}, {7, 1}, {7, 2}, {7, 3}, {7, 5}, {7, 6}, {8, 2}, {8, 3},
             {8, 5}, {8, 6}, {9, 0}, {9, 2}, {9, 3}, {9, 5}, {9, 6}, {9, 8}, {10, 1},
             {10, 2}, {10, 3}, {10, 5}, {10, 6}}}
));

TEST_P(GraphWriterToDimacsBinaryTest, test) {
    GraphIOParameters const& parameters = GetParam();
    NetworKit::Graph G(parameters.N, false, false);
    for (const auto &[u, v] : parameters.E) {
        G.addEdge(u, v);
    }
    std::string path = generate_filename("input");
    Koala::DimacsBinaryGraphWriter().write(G, path);
    EXPECT_TRUE(compare_files(path, parameters.G));
    remove(path.data());
}

INSTANTIATE_TEST_SUITE_P(
    test_small, GraphWriterToDimacsBinaryTest, testing::Values(
        GraphIOParameters{
            "input/example_1.col.b", 11,
            {{2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 2}, {4, 3}, {5, 0}, {5, 1}, {5, 3}, {5, 4},
             {6, 0}, {6, 1}, {6, 2}, {7, 1}, {7, 2}, {7, 3}, {7, 5}, {7, 6}, {8, 2}, {8, 3},
             {8, 5}, {8, 6}, {9, 0}, {9, 2}, {9, 3}, {9, 5}, {9, 6}, {9, 8}, {10, 1},
             {10, 2}, {10, 3}, {10, 5}, {10, 6}}}
));
