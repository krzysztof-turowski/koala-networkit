#include "shortest_path/Seidel.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

#include <networkit/components/ConnectedComponents.hpp>

namespace Koala {

static const int INF = std::numeric_limits<int>::max();

void SeidelAPSP::checkInput() const {
    if (graph->isDirected()) {
        throw std::invalid_argument("Graph must be undirected for Seidel's algorithm.");
    }
    if (graph->isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for Seidel's algorithm.");
    }
    // The APD recursion squares the graph until it becomes complete, which never
    // happens across components -> it would recurse forever on a disconnected
    // graph. Reject it up front instead of hanging.
    NetworKit::ConnectedComponents cc(*graph);
    cc.run();
    if (cc.numberOfComponents() > 1) {
        throw std::invalid_argument("Graph must be connected for Seidel's algorithm.");
    }
}

void SeidelAPSP::run() {
    const NetworKit::count N = graph->upperNodeIdBound();
    const int m = static_cast<int>(graph->numberOfNodes());

    std::vector<int> compact(N, -1);        // real id -> compact index
    std::vector<NetworKit::node> original;  // compact index -> real id
    original.reserve(m);
    for (const auto u : graph->nodeRange()) {
        compact[u] = static_cast<int>(original.size());
        original.push_back(u);
    }

    Matrix A = Matrix::Zero(m, m);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        A(compact[u], compact[v]) = 1;
        A(compact[v], compact[u]) = 1;
    });

    const Matrix D = APD(A);
    distances.assign(N, std::vector<int>(N, INF));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            distances[original[i]][original[j]] = D(i, j);
        }
    }

    computeDiameter();

    hasRun = true;
}

SeidelAPSP::Matrix SeidelAPSP::APD(const Matrix& A) const {
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
