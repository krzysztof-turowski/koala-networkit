// GoemansWilliamsonMaxCut.cpp
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>

#include <max_cut/GoemansWilliamsonMaxCut.hpp>

// temporary for debugging
#include <iostream>
#include <iomanip>

namespace Koala {

GoemansWilliamsonMaxCut::GoemansWilliamsonMaxCut(const std::vector<std::vector<int>>& graphInput)
    : graph(graphInput), n(graphInput.size()), maxCutValue(0) {
    srand(time(NULL));

    bestSet.resize(n, false);
}

void GoemansWilliamsonMaxCut::solve() {
    // Declarations
    struct blockmatrix C, X, Z;
    double *b = NULL, *y = NULL;
    double pobj, dobj;
    struct constraintmatrix *constraints = NULL;
    maxCutValue = 0;
    
    // Prepare constraints
    initializeSDP(C, b, constraints);

    // Calculate
    initsoln(n, 1, C, b, constraints, &X, &y, &Z);
    int ret = easy_sdp(n, 1, C, b, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);

    std::vector<double> r = randomUnitVector(n);
    
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += X.blocks[1].data.mat[i * n + j] * r[j];
        }
        bestSet[i] = (sum >= 0);
    }
    
    // Calculate cut
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (bestSet[i] != bestSet[j]) {
                maxCutValue += graph[i][j];
            }
        }
    }

    free_prob(n, 1, C, b, constraints, X, y, Z);
}

int GoemansWilliamsonMaxCut::getMaxCutValue() const {
    return maxCutValue;
}

const std::vector<bool>& GoemansWilliamsonMaxCut::getBestSet() const {
    return bestSet;
}

std::vector<double> GoemansWilliamsonMaxCut::randomUnitVector(int dim) {
    std::vector<double> r(dim);
    double norm = 0;
    for (int i = 0; i < dim; ++i) {
        r[i] = static_cast<double>(rand()) / RAND_MAX - 0.5;
        norm += r[i] * r[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < dim; ++i) {
        r[i] /= norm;
    }
    return r;
}

void GoemansWilliamsonMaxCut::initializeSDP(struct blockmatrix &C, double *&b, struct constraintmatrix *&constraints) {
    // Allocate the block matrix
    C.nblocks = 1;
    C.blocks = (struct blockrec *) malloc((C.nblocks + 1) * sizeof(struct blockrec));
    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = n;
    C.blocks[1].data.mat = (double *) malloc((n * n + 1) * sizeof(double));
    memset(C.blocks[1].data.mat, 0, (n * n + 1) * sizeof(double));

    // Set up the objective matrix C (weights)
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i != j && graph[i - 1][j - 1] != 0) {
                C.blocks[1].data.mat[ijtok(i, j, n)] = -graph[i - 1][j - 1];
                C.blocks[1].data.mat[ijtok(i, j, n)] = -graph[i - 1][j - 1];
            }
        }
    }

    // Linear constraints
    b = (double *) malloc(sizeof(double) * (n + 1));
    constraints = (struct constraintmatrix *) malloc((n + 1) * sizeof(struct constraintmatrix));
    for (int i = 1; i <= n; ++i) {
        b[i] = 1.0;
        constraints[i].blocks = (struct sparseblock *) malloc(sizeof(struct sparseblock));
        constraints[i].blocks->blocknum = 1;
        constraints[i].blocks->blocksize = n;   // should i add?
        constraints[i].blocks->constraintnum = i;
        constraints[i].blocks->entries = (double *) malloc(sizeof(double));
        constraints[i].blocks->iindices = (int *) malloc(sizeof(int));
        constraints[i].blocks->jindices = (int *) malloc(sizeof(int));
        constraints[i].blocks->numentries = 1;
        constraints[i].blocks->iindices[1] = i;
        constraints[i].blocks->jindices[1] = i;
        constraints[i].blocks->entries[1] = 1.0;
        constraints[i].blocks->next = NULL;
    }
}

} // namespace Koala
