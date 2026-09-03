#include "shortest_path/Seidel.hpp"

#include <stdexcept>

namespace Koala {

APDAlgorithm::APDAlgorithm(const NetworKit::Graph& G) : G(&G) {
    if (G.isDirected()) {
        throw std::invalid_argument("Graph must be undirected for APD algorithm.");
    }
    if (G.isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for APD algorithm.");
    }
}

void APDAlgorithm::run() {
    NetworKit::count n = G->upperNodeIdBound();

    Matrix A = Matrix::Zero(n, n);

    G->forEdges([&](NetworKit::node u, NetworKit::node v) {
        A(u, v) = 1;
        A(v, u) = 1;
    });

    distances = APD(A);
    hasRun = true;
}

const APDAlgorithm::Matrix& APDAlgorithm::getDistances() const {
    assureFinished();
    return distances;
}

int APDAlgorithm::getDistance(NetworKit::node u, NetworKit::node v) const {
    assureFinished();
    return distances(u, v);
}

APDAlgorithm::Matrix APDAlgorithm::APD(const Matrix& A) const {
    int n = A.rows();

    Matrix Z = A * A;

    Matrix B = Matrix::Zero(n, n);
    bool is_complete = true;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                if (A(i, j) == 1 || Z(i, j) > 0) {
                    B(i, j) = 1;
                } else {
                    is_complete = false;
                }
            }
        }
    }

    if (is_complete) {
        Matrix D = Matrix::Zero(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    D(i, j) = 2 * B(i, j) - A(i, j);
                }
            }
        }
        return D;
    }

    Matrix T = APD(B);
    Matrix X = T * A;

    Eigen::VectorXi degree = A.colwise().sum();

    Matrix D = Matrix::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                if (X(i, j) >= T(i, j) * degree(j)) {
                    D(i, j) = 2 * T(i, j);
                } else {
                    D(i, j) = 2 * T(i, j) - 1;
                }
            }
        }
    }

    return D;
}

}  // namespace Koala
