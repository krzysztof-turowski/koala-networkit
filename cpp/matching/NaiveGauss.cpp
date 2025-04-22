#include <Eigen/Core>
#include <Eigen/Dense>

#include <matching/utils.hpp>

#include <iostream>

using namespace std;
using namespace Eigen;

namespace Koala {
    bool eliminate(MatrixXd& A, int r, int c) {
        if (eq(A(r, c), 0)) return false;

        A -= A.col(c) * A.row(r) / A(r, c);
        return true;
    }

    vector<int> pivotElimination(MatrixXd& A, function<bool(int, int)> isCellAllowed) {
        int n = A.cols();

        vector<int> res(n);
        for (int c = 0; c < n; ++c) {
            for (int r = 0; r < n; ++r) {
                if (eq(A(r, c), 0) || !isCellAllowed(c, r))
                    continue;

                eliminate(A, r, c);
                res[c] = r;
                break;
            }
        }

        return res;
    }

    vector<int> simpleElimination(MatrixXd& A, int k) {
        vector<int> res;

        for (int i = 0; i < k; ++i) {
            if (eliminate(A, i, i)) {
                res.push_back(i);
            }
        };

        return res;
    }
}
