#pragma once

#include <networkit/graph/Graph.hpp>
#include <Eigen/Core>
#include <set>

typedef std::set<std::pair<int, int>> Matching;

namespace Koala {
    class GeneralGaussianMatching {
    public:
        GeneralGaussianMatching(const NetworKit::Graph& G);
        void run();
        Matching getMatching();

        // private:
        NetworKit::Graph G;
        Eigen::MatrixXd AG;
        Matching M;

        std::vector<int> oldIdx;
    };
}
