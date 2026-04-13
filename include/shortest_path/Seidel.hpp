#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <Eigen/Dense>

namespace Koala {

    class APDAlgorithm : public NetworKit::Algorithm {
    public:
        using Matrix = Eigen::MatrixXi;

        APDAlgorithm(const NetworKit::Graph& G);

        void run() override;

        const Matrix& getDistances() const;

        int getDistance(NetworKit::node u, NetworKit::node v) const;

    private:
        const NetworKit::Graph* G;
        Matrix distances;

        Matrix APD(const Matrix& A) const;
    };

} /* namespace koala */