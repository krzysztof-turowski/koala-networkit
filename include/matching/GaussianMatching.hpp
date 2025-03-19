#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <tuple>

typedef std::vector<std::pair<int, int>> Matching;

namespace Koala {
    class GaussianMatching : public NetworKit::Algorithm {
    public:
        GaussianMatching(NetworKit::Graph& graph, bool bipartite = false);

        Matching getMatching();

        void run();

    private:
        bool bipartite;
        NetworKit::Graph graph;

        Eigen::MatrixXd AG;
    };
}


// TODO: move to some utils

int generateRandom(int min, int max) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(min, max);

    return (int)dist6(rng);
}

constexpr double EPS = 1e-8;
bool eq0(double a) {
    return -EPS <= a && a <= EPS;
}