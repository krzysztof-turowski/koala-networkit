#include <gtest/gtest.h>

#include <list>

#include <independent_set/SimpleIndependentSet.hpp>

#include "helpers.hpp"

struct IndependentSetParameters {
    int N;
    std::list<std::pair<int, int>> E;
    int expectedSetSize;
};

template <typename Algorithm>
class SimpleGraphs : public testing::Test {
public:
    virtual void SetUp() {
    }

protected:
    void verify(IndependentSetParameters& parameters) {
        NetworKit::Graph G = build_graph(parameters.N, parameters.E);
        auto algorithm = Algorithm(G);
        algorithm.run();

        auto independentSet = algorithm.getIndependentSet();
        for (const auto &[u, v] : parameters.E) {
            EXPECT_FALSE(independentSet[u] && independentSet[v]);
        }

        int setSize = 0;
        for (const auto &[v, belongs] : independentSet) {
            setSize += (belongs);
        }
        EXPECT_EQ(setSize, parameters.expectedSetSize);
    }
};

/*
the following graphs are from:
https://mathworld.wolfram.com/IndependentSet.html
*/
TYPED_TEST_CASE_P(SimpleGraphs);

TYPED_TEST_P(SimpleGraphs, WheelGraphW_8) {
    IndependentSetParameters parameters =
    {8, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {7, 0}, {7, 1}, {7, 2}, {7, 3}, {7, 4}, {7, 5}, {7, 6},
        },
    3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, UtilityGraphK_3_3) {
    IndependentSetParameters parameters =
    {6, {
        {0, 3}, {0, 4}, {0, 5},
        {1, 3}, {1, 4}, {1, 5},
        {2, 3}, {2, 4}, {2, 5}
        },
    3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, PetersenGraph) {
    IndependentSetParameters parameters =
    {10, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
        {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9},
        {5, 7}, {5, 8}, {6, 8}, {6, 9}, {7, 9}
        },
    4};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, FruchtGraph) {
    IndependentSetParameters parameters =
    {12, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 0},
        {0, 7}, {1, 8}, {2, 8}, {3, 9}, {4, 9}, {5, 10}, {6, 10},
        {7, 8}, {11, 7}, {11, 9}, {11, 10}
    },
    5};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, TroubleGraph_1) {
    IndependentSetParameters parameters =
    {7, {
        {0,2}, {0,3}, {1,3}, {0,4}, {1,4}, {3,4}, {0,5}, {1,5}, {2,5}, {0,6}, {1,6}, {2,6},
    },
    3};
    this->verify(parameters);
}

TYPED_TEST_P(SimpleGraphs, TwoK5) {
    IndependentSetParameters parameters =
    {10, {
        {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4},
        {5,6}, {5,7}, {5,8}, {5,9}, {6,7}, {6,8}, {6,9}, {7,8}, {7,9}, {8,9},
    },
    2};
    this->verify(parameters);
}
REGISTER_TYPED_TEST_CASE_P(SimpleGraphs, WheelGraphW_8, UtilityGraphK_3_3, PetersenGraph, FruchtGraph, TroubleGraph_1, TwoK5);

typedef testing::Types<
    Koala::BruteForceIndependentSet,
    Koala::Mis1IndependentSet,
    Koala::Mis2IndependentSet,
    Koala::Mis3IndependentSet,
    Koala::Mis4IndependentSet,
    Koala::Mis5IndependentSet
    >Algorithms;
INSTANTIATE_TYPED_TEST_CASE_P(IndependentSet, SimpleGraphs, Algorithms);


class CompareAlgorithmResults : public testing::Test {
public:
    virtual void SetUp() {
    }

protected:

};

template <typename Algorithm>
std::optional<int> runAndValidate(NetworKit::Graph G, std::list<std::pair<int, int>>& edges) {
    auto algorithm = Algorithm(G);
    algorithm.run();

    auto independentSet = algorithm.getIndependentSet();
    for (const auto &[u, v] : edges) {
        bool neighbors = independentSet[u] && independentSet[v];
        EXPECT_FALSE(neighbors);
        if (neighbors) {
            return std::nullopt;
        }
    }

    int setSize = 0;
    for (const auto &[v, belongs] : independentSet) {
        setSize += (belongs);
    }
    return setSize;
}

class AdjacencyMatrix { // assumes that graph size at most 11
public:
    AdjacencyMatrix(int numberOfVertices) : numberOfVertices(numberOfVertices) {
        numberOfEdges = numberOfVertices * (numberOfVertices - 1) / 2;
        assert(numberOfEdges <= 64);
        biggestRepresentation = (1ull << numberOfEdges) - 1;
    };
    bool hasNext() {
        return representation <= biggestRepresentation;
    };

    void goNext() {
        ++representation;
    }

    std::list<std::pair<int, int>> getNiceEdges() {
        std::list<std::pair<int, int>> list;
        for (int i = 0; i < numberOfEdges; ++i) {
            if (representation & (1 << i)) {
                list.push_back(edges[i]);
            }
        }
        return list;
    }

private:
    std::pair<int, int> edges[55] {
        {0, 1},
        {0, 2}, {1, 2},
        {0, 3}, {1, 3}, {2, 3},
        {0, 4}, {1, 4}, {2, 4}, {3, 4},
        {0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5},
        {0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6},
        {0, 7}, {1, 7}, {2, 7}, {3, 7}, {4, 7}, {5, 7}, {6, 7},
        {0, 8}, {1, 8}, {2, 8}, {3, 8}, {4, 8}, {5, 8}, {6, 8}, {7, 8},
        {0, 9}, {1, 9}, {2, 9}, {3, 9}, {4, 9}, {5, 9}, {6, 9}, {7, 9}, {8, 9},
        {0, 10}, {1, 10}, {2, 10}, {3, 10}, {4, 10}, {5, 10}, {6, 10}, {7, 10}, {8, 10}, {9, 10},
    };

    int numberOfVertices;
    int numberOfEdges;
    unsigned long long representation = 0;
    unsigned long long biggestRepresentation;
};

TEST(CompareAlgorithmResults, test) {
    constexpr int maximumGraphSize = 6;
    for (int numberOfVertices = 1; numberOfVertices <= maximumGraphSize; ++numberOfVertices) {
        AdjacencyMatrix adj(numberOfVertices);
        while(true) {
            std::list<std::pair<int, int>> edges = adj.getNiceEdges();
            NetworKit::Graph G = build_graph(numberOfVertices, edges);

            std::vector<std::optional<int>> algorithmSetSizes {
                //runAndValidate<Koala::BruteForceIndependentSet>(G, edges),
                //runAndValidate<Koala::Mis1IndependentSet>(G, edges),
                runAndValidate<Koala::Mis2IndependentSet>(G, edges),
                runAndValidate<Koala::Mis3IndependentSet>(G, edges),
                //runAndValidate<Koala::Mis4IndependentSet>(G, edges),
                //runAndValidate<Koala::Mis5IndependentSet>(G, edges),
            };

            int bestAchievedSize = 0;
            for (auto size : algorithmSetSizes) {
                if (size.has_value()) {
                    bestAchievedSize = std::max(bestAchievedSize, *size);
                }
                else {
                    return;
                }
            }

            for (int i = 0; i < algorithmSetSizes.size(); ++i) {
                EXPECT_TRUE(algorithmSetSizes[i].has_value());
                if (algorithmSetSizes[i].has_value()) {
                    EXPECT_EQ(algorithmSetSizes[i], bestAchievedSize);
                    if (algorithmSetSizes[i] != bestAchievedSize) {
                        for (auto e : edges) {
                            std::cout << "(" << e.first << "," << e.second << ") ";
                        } 
                        std::cout << std::endl;  
                        return;
                    }
                }
            }

            if (!adj.hasNext()) {

                break;
            }
            adj.goNext();
        }
    }
}
