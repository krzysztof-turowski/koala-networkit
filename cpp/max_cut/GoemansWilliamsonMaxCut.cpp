/*
 * GoemansWilliamsonMaxCut.cpp
 *
 * Solution for the Max-Cut problem using the Goemans-Williamson algorithm.
 * Created on: 22.04.2024
 * Author: Michał Miziołek
 */

extern "C" {
#include <declarations.h>
}

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <random>
#include <stdexcept>
#include <vector>

#include <max_cut/GoemansWilliamsonMaxCut.hpp>

namespace Koala {
namespace {

NetworKit::edgeweight flipGain(
        const NetworKit::Graph &graph, const std::vector<bool> &cut, NetworKit::node u) {
    NetworKit::edgeweight gain = 0;
    graph.forNeighborsOf(u, [&](NetworKit::node v, NetworKit::edgeweight weight) {
        gain += cut[u] == cut[v] ? weight : -weight;
    });
    return gain;
}

void improveCut(const NetworKit::Graph &graph, std::vector<bool> &cut) {
    bool improved;
    do {
        improved = false;
        for (auto u : graph.nodeRange()) {
            if (flipGain(graph, cut, u) > 0) {
                cut[u] = !cut[u];
                improved = true;
            }
        }
    } while (improved);
}

}  // namespace

void GoemansWilliamsonMaxCut::run() {
    // Declarations
    int n = graph->numberOfNodes();
    maxCutSet.assign(n, false);
    maxCutValue = 0;
    if (n == 0) {
        hasRun = true;
        return;
    }

    struct blockmatrix C, X, Z;
    double *b = NULL, *y = NULL;
    double pobj, dobj;
    struct constraintmatrix *constraints = NULL;

    // Prepare constraints
    initializeSDP(C, b, constraints);

    // Calculate
    initsoln(n, 1, C, b, constraints, &X, &y, &Z);
    easy_sdp(n, 1, C, b, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);

    Eigen::MatrixXd covariance(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            covariance(i, j) = X.blocks[1].data.mat[ijtok(i + 1, j + 1, n)];
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
    if (solver.info() != Eigen::Success) {
        free_prob(n, 1, C, b, constraints, X, y, Z);
        throw std::runtime_error("Failed to factor the MaxCut SDP solution");
    }
    Eigen::VectorXd eigenvalues = solver.eigenvalues().cwiseMax(0).cwiseSqrt();
    Eigen::MatrixXd vectors = solver.eigenvectors() * eigenvalues.asDiagonal();

    const int rounds = std::clamp(4 * n, 32, 128);
    for (int round = 0; round < rounds; ++round) {
        std::vector<double> random_vector = randomUnitVector(n);
        Eigen::Map<Eigen::VectorXd> hyperplane(random_vector.data(), n);
        std::vector<bool> cut(n);
        for (int i = 0; i < n; ++i) {
            cut[i] = vectors.row(i).dot(hyperplane) >= 0;
        }
        improveCut(*graph, cut);

        double cut_value = calculateCutValue(cut);
        if (cut_value > maxCutValue) {
            maxCutValue = cut_value;
            maxCutSet = std::move(cut);
        }
    }

    free_prob(n, 1, C, b, constraints, X, y, Z);
    hasRun = true;
}

std::vector<double> GoemansWilliamsonMaxCut::randomUnitVector(int dim) {
    std::vector<double> r(dim);
    double norm = 0;
    std::mt19937 generator(std::random_device {}());
    std::normal_distribution<double> distribution;

    for (int i = 0; i < dim; ++i) {
        r[i] = distribution(generator);
        norm += r[i] * r[i];
    }

    norm = std::sqrt(norm);
    for (int i = 0; i < dim; ++i) {
        r[i] /= norm;
    }

    return r;
}

void GoemansWilliamsonMaxCut::initializeSDP(struct blockmatrix &C, double *&b,
                struct constraintmatrix *&constraints) {
    int n = graph->numberOfNodes();

    // Allocate the block matrix
    C.nblocks = 1;
    C.blocks = reinterpret_cast<struct blockrec *>
                    (malloc((C.nblocks + 1) * sizeof(struct blockrec)));
    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = n;
    C.blocks[1].data.mat = reinterpret_cast<double *>(malloc((n * n + 1) * sizeof(double)));
    memset(C.blocks[1].data.mat, 0, (n * n + 1) * sizeof(double));

    // Set up the objective matrix C (weights)
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        C.blocks[1].data.mat[ijtok(u + 1, v + 1, n)] = -w;
        C.blocks[1].data.mat[ijtok(v + 1, u + 1, n)] = -w;
    });

    // Linear constraints
    b = reinterpret_cast<double *>(malloc(sizeof(double) * (n + 1)));
    constraints = reinterpret_cast<struct constraintmatrix *>
                    (malloc((n + 1) * sizeof(struct constraintmatrix)));
    for (int i = 1; i <= n; ++i) {
        b[i] = 1.0;
        constraints[i].blocks = reinterpret_cast<struct sparseblock *>
                    (malloc(sizeof(struct sparseblock)));
        constraints[i].blocks->blocknum = 1;
        constraints[i].blocks->blocksize = n;
        constraints[i].blocks->constraintnum = i;
        constraints[i].blocks->entries = reinterpret_cast<double*>(malloc(2 * sizeof(double)));
        constraints[i].blocks->iindices = reinterpret_cast<int*>(malloc(2 * sizeof(int)));
        constraints[i].blocks->jindices = reinterpret_cast<int*>(malloc(2 * sizeof(int)));
        constraints[i].blocks->numentries = 1;
        constraints[i].blocks->iindices[1] = i;
        constraints[i].blocks->jindices[1] = i;
        constraints[i].blocks->entries[1] = 1.0;
        constraints[i].blocks->next = NULL;
    }
}

}  // namespace Koala
