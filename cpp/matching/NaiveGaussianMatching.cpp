// #include <matching/GaussianMatching.hpp>

// #include<iostream>

// using namespace Eigen;
// using namespace NetworKit;

// int generateRandom(int min, int max);
// bool eq(double a, double b);

// namespace Koala {

//     GaussianMatching::GaussianMatching(Graph& G) : GaussianMatching(G) {}

//     Matching GaussianMatching::getMatching() {
//         return M;
//     }

//     void NaiveGaussianMatching::run() {
//         int n = G.numberOfNodes();
//         M = {};

//         double det = AG.determinant();
//         if (eq(det, 0)) return;

//         MatrixXd B = AG.inverse();
//         for (int c = 0; c < n; ++c) {
//             int r;
//             for (r = 0; r < n; ++r) {
//                 if (eq(B(r, c), 0) || eq(AG(c, r), 0)) continue;

//                 bool eliminated = true;
//                 for (int i = 0; i < n; ++i) {
//                     if (i == r) continue;
//                     if (eq(B(r, i), 0)) {
//                         eliminated = false;
//                         break;
//                     }
//                 }
//                 if (!eliminated) break;
//             }
//             if (r == n) return;

//             B -= B.col(c) * B.row(r) / B(r, c);
//             // if (!bipartite) {
//             //     B -= B.col(r) * B.row(c) / B(c, r);
//             // }

//             M.push_back({ c + 1, r + 1 });
//         }
//     }
// }


// // TODO: move to some utils

// int generateRandom(int min, int max) {
//     std::random_device dev;
//     std::mt19937 rng(dev());
//     std::uniform_int_distribution<std::mt19937::result_type> dist6(min, max);

//     return (int)dist6(rng);
// }

// constexpr double EPS = 1e-8;
// bool eq(double a, double b) {
//     return fabs(a - b) <= EPS;
// }