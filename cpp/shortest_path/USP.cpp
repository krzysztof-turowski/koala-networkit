#include "shortest_path/USP.hpp"

#include <cmath>
#include <stdexcept>

namespace Koala {

// Element-wise (.array()) boolean helpers for Eigen Matrix
namespace {

USPAlgorithm::Matrix elemAnd(const USPAlgorithm::Matrix& A, const USPAlgorithm::Matrix& B) {
    return (A.array() * B.array()).matrix();
}

USPAlgorithm::Matrix elemOr(const USPAlgorithm::Matrix& A, const USPAlgorithm::Matrix& B) {
    return (A.array() + B.array() - A.array() * B.array()).matrix();
}

USPAlgorithm::Matrix elemNot(const USPAlgorithm::Matrix& A) {
    return (1 - A.array()).matrix();
}

}

USPAlgorithm::USPAlgorithm(const NetworKit::Graph& G) : G(&G) {
    if (G.isDirected()) {
        throw std::invalid_argument("Graph must be undirected for USP algorithm.");
    }
    if (G.isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for USP algorithm.");
    }
}

USPAlgorithm::Matrix USPAlgorithm::boolMul(const Matrix& A, const Matrix& B) {
    Matrix C = A * B;
    return (C.array() > 0).cast<int>();
}

void USPAlgorithm::run() {
    const int n = static_cast<int>(G->upperNodeIdBound());

    int l = 0;
    while ((1 << l) < n) ++l;

    Matrix A = Matrix::Zero(n, n);
    G->forEdges([&](NetworKit::node u, NetworKit::node v) {
        A(u, v) = 1;
        A(v, u) = 1;
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

    std::vector<Matrix> Ck(l + 1, Matrix::Zero(n, n));
    std::vector<Matrix> Dk(l + 1, Matrix::Zero(n, n));
    std::vector<Matrix> Pk(l + 1, Matrix::Zero(n, n));
    std::vector<Matrix> Qk(l + 1, Matrix::Zero(n, n));

    Ck[l] = Matrix::Ones(n, n);
    Dk[l] = Matrix::Ones(n, n);
    Pk[l] = Matrix::Identity(n, n);

    for (int k = l - 1; k >= 0; --k) {
        Matrix notC = elemNot(Ck[k + 1]);
        Matrix PA = boolMul(Pk[k + 1], Ak[k]);
        Matrix QA = boolMul(Qk[k + 1], Ak[k]);

        Ck[k] = elemOr(elemAnd(PA, Ck[k + 1]), elemAnd(QA, notC));

        Matrix PB = boolMul(Pk[k + 1], Bk[k]);
        Matrix QB = boolMul(Qk[k + 1], Bk[k]);

        Dk[k] = elemOr(elemAnd(PB, Ck[k + 1]), elemAnd(QB, notC));

        Pk[k] = elemOr(Pk[k + 1], Qk[k + 1]);
        Qk[k] = elemAnd(Dk[k], elemNot(Ck[k]));
    }

    distances = Matrix::Zero(n, n);
    for (int k = 0; k <= l; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (Ck[k](i, j) == 0) {
                    distances(i, j) += (1 << k);
                }
            }
        }
    }

    hasRun = true;
}

const USPAlgorithm::Matrix& USPAlgorithm::getDistances() const {
    assureFinished();
    return distances;
}

int USPAlgorithm::getDistance(NetworKit::node u, NetworKit::node v) const {
    assureFinished();
    return distances(u, v);
}

}
