#include <matching/GaussianMatching.hpp>

using namespace Eigen;
using namespace NetworKit;

namespace Koala {

    GaussianMatching::GaussianMatching(Graph& graph, bool bipartite=false):
        bipartite(bipartite) {
        int n = graph.numberOfNodes();
        AG = ArrayXXd::Zero(n, n);

        graph.forEdges([this](
            node u, node v,
            edgeweight weight, edgeid id) {
                int Xuv = generateRandom(1, 1000);
                int minV = std::min(u, v), maxV = std::max(u, v);
                AG(minV, maxV) = Xuv;
                AG(maxV, minV) = -Xuv;
            });
    }

    Matching GaussianMatching::getMatching() {
        int n = graph.numberOfNodes();

        MatrixXd B = AG.inverse();
        Matching M = {};

        for (int c = 0; c < n; ++c) {
            int r;
            for (r = 0; r < n; ++r) {
                if (eq0(B(r, c)) || eq0(AG(c,r))) continue;

                bool eliminated = true;
                for (int i = 0; i < n; ++i) {
                    if (i == r) continue;
                    if (B(r, i) != 0) {
                        eliminated = false;
                        break;
                    }
                }
                if (!eliminated) break;
            }

            B -= B.col(c) * B.row(r) / B(c,r);
            if (!bipartite) {
                B -= B.col(r) * B.row(c) / B(r,c);
            }    

            M.push_back({c,r});
        }

        return M;
    }

    void GaussianMatching::run() {
        
    }
}
