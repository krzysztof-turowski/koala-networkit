#include "shortest_path/ShoshanZwick.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/components/ConnectedComponents.hpp>

namespace Koala {

static const int INF = std::numeric_limits<int>::max();

// Element-wise (.array()) boolean helpers for Eigen Matrix
namespace {

ShoshanZwickAPSP::Matrix elemAnd(const ShoshanZwickAPSP::Matrix& A, const ShoshanZwickAPSP::Matrix& B) {
    return (A.array() * B.array()).matrix();
}

ShoshanZwickAPSP::Matrix elemOr(const ShoshanZwickAPSP::Matrix& A, const ShoshanZwickAPSP::Matrix& B) {
    return (A.array() + B.array() - A.array() * B.array()).matrix();
}

ShoshanZwickAPSP::Matrix elemNot(const ShoshanZwickAPSP::Matrix& A) {
    return (1 - A.array()).matrix();
}

}

void ShoshanZwickAPSP::checkInput() const {
    if (graph->isDirected()) {
        throw std::invalid_argument("Graph must be undirected for Shoshan-Zwick's algorithm.");
    }
    if (graph->isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for Shoshan-Zwick's algorithm.");
    }
    // On a disconnected graph the bit-recovery never marks cross-component pairs
    // as unreachable, so it returns wrong finite distances instead of infinities.
    NetworKit::ConnectedComponents cc(*graph);
    cc.run();
    if (cc.numberOfComponents() > 1) {
        throw std::invalid_argument("Graph must be connected for Shoshan-Zwick's algorithm.");
    }
}

ShoshanZwickAPSP::Matrix ShoshanZwickAPSP::boolMul(const Matrix& A, const Matrix& B) {
    Matrix C = A * B;
    return (C.array() > 0).cast<int>();
}

void ShoshanZwickAPSP::run() {
    const NetworKit::count N = graph->upperNodeIdBound();
    const int n = static_cast<int>(graph->numberOfNodes());

    std::vector<int> compact(N, -1);        // real id -> compact index
    std::vector<NetworKit::node> original;  // compact index -> real id
    original.reserve(n);
    for (const auto u : graph->nodeRange()) {
        compact[u] = static_cast<int>(original.size());
        original.push_back(u);
    }

    int l = 0;
    while ((1 << l) < n) ++l;

    Matrix A = Matrix::Zero(n, n);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        A(compact[u], compact[v]) = 1;
        A(compact[v], compact[u]) = 1;
    });

    std::vector<Matrix> Ak(l + 1);
    std::vector<Matrix> Bk(l + 1);

    Ak[0] = Matrix::Identity(n, n);
    Bk[0] = A;
    for (int i = 0; i < n; ++i) {
        Bk[0](i, i) = 1;
    }

    for (int k = 1; k <= l; ++k) {
        Ak[k] = boolMul(Ak[k - 1], Bk[k - 1]);
        Bk[k] = boolMul(Bk[k - 1], Bk[k - 1]);
    }

    Matrix C = Matrix::Ones(n, n);
    Matrix P = Matrix::Identity(n, n);
    Matrix Q = Matrix::Zero(n, n);

    distances.assign(N, std::vector<int>(N, INF));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            distances[original[i]][original[j]] = 0;
        }
    }

    for (int k = l - 1; k >= 0; --k) {
        Matrix notC = elemNot(C);
        Matrix PA = boolMul(P, Ak[k]);
        Matrix QA = boolMul(Q, Ak[k]);
        Matrix Cnext = elemOr(elemAnd(PA, C), elemAnd(QA, notC));

        Matrix PB = boolMul(P, Bk[k]);
        Matrix QB = boolMul(Q, Bk[k]);
        Matrix Dnext = elemOr(elemAnd(PB, C), elemAnd(QB, notC));

        Matrix Pnext = elemOr(P, Q);
        Matrix Qnext = elemAnd(Dnext, elemNot(Cnext));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (Cnext(i, j) == 0) {
                    distances[original[i]][original[j]] += (1 << k);
                }
            }
        }

        C = std::move(Cnext);
        P = std::move(Pnext);
        Q = std::move(Qnext);
    }

    computeDiameter();

    hasRun = true;
}

}  // namespace Koala
