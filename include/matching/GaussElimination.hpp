#pragma once

#include <Eigen/Core>

namespace Koala {
    void eliminate(Eigen::MatrixXd& A, int r, int c);

    std::vector<int> pivotElimination(Eigen::MatrixXd& A, std::function<bool(int, int)> isCellAllowed);

    std::vector<int> simpleElimination(Eigen::MatrixXd& A, int k);
}