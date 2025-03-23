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
        Matching M;
    };
}
