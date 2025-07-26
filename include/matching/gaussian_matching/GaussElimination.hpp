#pragma once

#include <Eigen/Core>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <matching/gaussian_matching/utils.hpp>

namespace Koala {
    void eliminate(MatZp& A, int r, int c);

    std::vector<int> pivotElimination(MatZp& A, std::function<bool(int, int)> isCellAllowed);

    std::vector<int> simpleElimination(MatZp& A, int k);
}