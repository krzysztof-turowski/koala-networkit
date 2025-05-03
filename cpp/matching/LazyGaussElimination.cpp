#include <Eigen/Core>
#include <Eigen/Dense>

#include <matching/utils.hpp>

#include <iostream>

using namespace std;
using namespace Eigen;

namespace Koala {
    typedef tuple<VectorXd, VectorXd, double> lazyElimination;

    static pair<int, int> get2div(int n) {
        // return m such that largest m: 2^m|n
        int m = 0, m2 = 1;
        while (n % (m2 * 2) == 0)
            m = m + 1, m2 *= 2;
        return { m, m2 };
    }

    static MatrixXd update(
        int r1, int r2, int c1, int c2,
        vector<lazyElimination> acc
    ) {
        int k = acc.size();
        MatrixXd U(r2 - r1 + 1, k);
        MatrixXd V(k, c2 - c1 + 1);
        for (int i = 0; i < k; ++i) {
            auto [u, v, a] = acc[i];
            U.col(i) = u.segment(r1, r2 - r1 + 1) / a;
            V.row(i) = v.segment(c1, c2 - c1 + 1);
        }
        return U * V;
    }

    vector<int> pivotElimination(MatrixXd& A, function<bool(int, int)> isCellAllowed) {
        int n = A.cols();

        vector<tuple<VectorXd, double>>lazy(n);
        vector<int> order(n);
        for (int c = 0; c < n; ++c) {
            int r;
            for (r = 0; r < n; ++r) {
                if (!eq(A(r, c), 0) && isCellAllowed(c, r))
                    break;
            }
            assert(r != n);
            order[c] = r;

            auto [j, j2] = get2div(c + 1);

            auto a = A(r, c);
            lazy[c] = { A.col(c), a };

            vector<tuple<VectorXd, VectorXd, double>> superLazy(j2);
            for (int i = 0; i < j2; ++i) {
                auto [l, l2] = get2div(i + 1);

                auto v = A.row(order[c - j2 + 1 + i]);
                auto [u, a] = lazy[c - j2 + 1 + i];

                auto from = superLazy.begin() + max(0, i - l2 + 1);
                auto to = superLazy.begin() + i; // todo
                auto acc = vector<lazyElimination>(from, to);

                int r1 = 0, r2 = n - 1; // todo
                int c1 = c + 1, c2 = min(c + j2, n - 1);

                if (r1 <= r2 && c1 <= c2) {
                    auto B = update(0, n - 1, c1, c2, acc);
                    A.block(r1, c1, r2 - r1 + 1, c2 - c1 + 1) -= B;
                }

                auto v1 = A.row(order[c - j2 + 1 + i]);
                auto [u1, a1] = lazy[c - j2 + 1 + i];
                superLazy[i] = { u1,v1,a1 };
            }

            A.col(c).setZero();
            A(r, c) = a;
        }
        return order;
    }


    vector<int> simpleElimination(MatrixXd& A, int k) {
        int n = A.cols();
        vector<int> res;

        vector<tuple<VectorXd, VectorXd, double>>lazy(k);
        for (int i = 0; i < k; ++i) {
            if (eq(A(i, i), 0)) {
                lazy[i] = { VectorXd::Zero(A.cols()), VectorXd::Zero(A.rows()), 1 };
            } else {
                lazy[i] = { A.col(i), A.row(i), A(i,i) };
                res.push_back(i);
            }
            auto a = A(i, i);

            auto [j, j2] = get2div(i + 1);

            auto from = lazy.begin() + max(0, i - j2 + 1);
            auto to = lazy.begin() + i + 1;
            auto acc = vector<lazyElimination>(from, to);

            int r1 = i + 1, r2 = min(i + j2, n - 1);
            int c1 = i + 1, c2 = n - 1;
            if (r1 <= r2 && c1 <= c2) {
                auto B = update(r1, r2, c1, c2, acc);
                A.block(r1, c1, r2 - r1 + 1, c2 - c1 + 1) -= B;
            }

            r1 = i + j2 + 1, r2 = n - 1;
            c1 = i + 1, c2 = min(i + j2, n - 1);
            if (r1 <= r2 && c1 <= c2) {
                auto B = update(r1, r2, c1, c2, acc);
                A.block(r1, c1, r2 - r1 + 1, c2 - c1 + 1) -= B;
            }

            A.col(i).setZero();
            A.row(i).setZero();
            A(i, i) = a;
        }

        return res;
    }
}
