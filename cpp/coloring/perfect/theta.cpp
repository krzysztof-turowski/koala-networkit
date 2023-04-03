#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" {
#include <declarations.h>
}

double theta(int n, int m, int *from, int *to) {
    int i, j, start, finish, temp;
    double pobj, dobj, *y, *a;
    struct blockmatrix C, X, Z;
    struct constraintmatrix *constraints;

    C.nblocks = 1;
    C.blocks = (blockrec*) malloc(2 * sizeof(struct blockrec));
    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = n;
    C.blocks[1].data.mat = (double*) malloc(n * n * sizeof(double));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            C.blocks[1].data.mat[ijtok(i, j, n)] = 1.0;
        }
    }
    constraints = (struct constraintmatrix*) malloc((m + 2) * sizeof(struct constraintmatrix));
    a = (double*) malloc((m + 2)* sizeof(double));
    a[1] = 1.0;
    constraints[1].blocks = (struct sparseblock*) malloc(sizeof(struct sparseblock));
    constraints[1].blocks->blocknum = constraints[1].blocks->constraintnum = 1;
    constraints[1].blocks->numentries = constraints[1].blocks->blocksize = n;
    constraints[1].blocks->next = constraints[1].blocks->nextbyblock = NULL;
    constraints[1].blocks->entries = (double*) malloc((n + 1) * sizeof(double));
    constraints[1].blocks->iindices = (int*) malloc((n + 1) * sizeof(int));
    constraints[1].blocks->jindices = (int*) malloc((n + 1) * sizeof(int));
    for (i = 1; i <= n; i++) {
        constraints[1].blocks->entries[i] = 1.0;
        constraints[1].blocks->iindices[i] = constraints[1].blocks->jindices[i] = i;
    }
    for (i = 2; i <= m + 1; i++) {
        a[i] = 0.0;
        constraints[i].blocks = (struct sparseblock*) malloc(sizeof(struct sparseblock));
        constraints[i].blocks->blocknum = constraints[i].blocks->numentries = 1;
        constraints[i].blocks->blocksize = n;
        constraints[i].blocks->constraintnum = i;
        constraints[i].blocks->next = constraints[i].blocks->nextbyblock = NULL;
        constraints[i].blocks->entries = (double*) malloc((2)* sizeof(double));
        constraints[i].blocks->iindices = (int*) malloc((2)* sizeof(int));
        constraints[i].blocks->jindices = (int*) malloc((2)* sizeof(int));
        constraints[i].blocks->entries[1] = 1.0;
        start = from[i - 2], finish = to[i - 2];
        if (start > finish) {
            std::swap(start, finish);
        }
        constraints[i].blocks->iindices[1] = start;
        constraints[i].blocks->jindices[1] = finish;
    }
    initsoln(n, m + 1, C, a, constraints, &X, &y, &Z);
    int ret = easy_sdp(n, m + 1, C, a, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);
    free_prob(n, m + 1, C, a, constraints, X, y, Z);
    return (dobj + pobj) / 2;
}
