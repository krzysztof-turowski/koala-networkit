/*
 * GoemansWilliamsonMaxCut.cpp
 *
 * Solution for the Max-Cut problem using the Goemans-Williamson algorithm.
 * Created on: 22.04.2024
 * Author: Michał Miziołek
 */

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <random>

extern "C" {
#include <declarations.h>
}

#include <max_cut/GoemansWilliamsonMaxCut.hpp>

namespace Koala {

void GoemansWilliamsonMaxCut::run() {
    // Declarations
    int n = graph->numberOfNodes();
    maxCutSet.assign(n, false);
    struct blockmatrix C, X, Z;
    double *b = NULL, *y = NULL;
    double pobj, dobj;
    struct constraintmatrix *constraints = NULL;
    maxCutValue = 0;

    // Prepare constraints
    initializeSDP(C, b, constraints);

    // Calculate
    initsoln(n, 1, C, b, constraints, &X, &y, &Z);
    easy_sdp(n, 1, C, b, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);

    std::vector<double> r = randomUnitVector(n);

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += X.blocks[1].data.mat[i * n + j] * r[j];
        }
        maxCutSet[i] = (sum >= 0);
    }

    // Calculate cut
    maxCutValue = calculateCutValue(maxCutSet);

    free_prob(n, 1, C, b, constraints, X, y, Z);
}

std::vector<double> GoemansWilliamsonMaxCut::randomUnitVector(int dim) {
    std::vector<double> r(dim);
    double norm = 0;
    std::mt19937 generator(std::random_device {}());
    std::uniform_real_distribution<double> distribution(-0.5, 0.5);

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
        constraints[i].blocks->entries = reinterpret_cast<double *>(malloc(sizeof(double)));
        constraints[i].blocks->iindices = reinterpret_cast<int *>(malloc(sizeof(int)));
        constraints[i].blocks->jindices = reinterpret_cast<int *>(malloc(sizeof(int)));
        constraints[i].blocks->numentries = 1;
        constraints[i].blocks->iindices[1] = i;
        constraints[i].blocks->jindices[1] = i;
        constraints[i].blocks->entries[1] = 1.0;
        constraints[i].blocks->next = NULL;
    }
}

}  // namespace Koala
