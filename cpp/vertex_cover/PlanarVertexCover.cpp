#include "vertex_cover/PlanarVertexCover.hpp"
#include <bit>
#include <queue>
#include <vector>
#include "networkit/Globals.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "techniques/separator/MISP.hpp"
#include "vertex_cover/VertexCover.hpp"

namespace Koala {
PlanarSeparatorVertexCover::PlanarSeparatorVertexCover(const NetworKit::Graph &G)
    : VertexCover(G) {}

// It is worth noting that in contrary to prepare algorithm described in **Bar-Yehuda, Reuven, and
// Shimon Even. "On approximating a vertex cover for planar graphs." Proceedings of the fourteenth
// annual ACM symposium on Theory of computing. 1982.** which is linear, this algorithm is O(n log
// n). In the end it does not affect the time complexity of the whole vertex cover algorithm, as
// usage of MISP algorithm is O(n log n) itself.
void PlanarSeparatorVertexCover::prepare(const NetworKit::Graph &G, std::vector<bool> &U,
                                         std::vector<bool> &VC, int n) {
    bool stop = false;

    while (!stop) {
        std::vector<int> deg_residual(G.upperNodeIdBound());
        G.forNodes([&](NetworKit::node v) {
            if (U[v])
                return;
            int deg = 0;
            G.forNeighborsOf(v, [&](NetworKit::node t) {
                if (!U[t])
                    deg++;
            });
            deg_residual[v] = deg;
        });

        std::queue<NetworKit::node> Q;
        G.forNodes([&](NetworKit::node v) {
            if (!U[v] && deg_residual[v] < 2)
                Q.push(v);
        });

        while (!Q.empty()) {
            auto v = Q.front();
            Q.pop();

            if (U[v])
                continue;

            if (deg_residual[v] == 0) {
                U[v] = true;
            }
            if (deg_residual[v] == 1) {
                NetworKit::node u = NetworKit::none;
                G.forNeighborsOf(v, [&](NetworKit::node t) {
                    if (u == NetworKit::none && !U[t])
                        u = t;
                });
                U[v] = true;
                U[u] = true;
                VC[u] = true;
                G.forNeighborsOf(u, [&](NetworKit::node t) {
                    if (!U[t] && --deg_residual[t] < 2)
                        Q.push(t);
                });
            }
        }

        // here deg(v) >= 2 for each v
        std::vector<bool> visited(G.upperNodeIdBound(), true);
        G.forNodes([&](NetworKit::node t) {
            if (!U[t])
                visited[t] = false;
        });

        for (int i = 0; i < static_cast<int>(visited.size()); i++) {
            if (!visited[i]) {
                std::vector<std::pair<NetworKit::node, int>> seen_deg_pair;
                Q.push(i);
                visited[i] = true;

                while (!Q.empty()) {
                    auto v = Q.front();
                    Q.pop();

                    seen_deg_pair.push_back({v, deg_residual[v]});

                    G.forNeighborsOf(v, [&](NetworKit::node t) {
                        if (!(visited[t] || U[t])) {
                            Q.push(t);
                            visited[t] = true;
                        }
                    });
                }

                bool all_degs_two = true;
                for (auto [v, deg] : seen_deg_pair) {
                    if (deg != 2) {
                        all_degs_two = false;
                        break;
                    }
                }

                if (all_degs_two) {
                    size_t len = seen_deg_pair.size();
                    auto cur = seen_deg_pair.front().first;
                    auto prev = NetworKit::none;

                    for (size_t k = 0; k < len; k++) {
                        if (k % 2 == 0)
                            VC[cur] = true;

                        NetworKit::node next = NetworKit::none;
                        for (auto v : G.neighborRange(cur)) {
                            if (v != prev && !U[v]) {
                                next = v;
                                break;
                            }
                        }

                        U[cur] = true;
                        prev = cur;
                        cur = next;
                    }
                }
            }
        }

        // bipartite finding
        auto bipGraph = bipartite(G, U, deg_residual);

        size_t UTrueSize = 0;
        for (auto v : U) {
            if (v)
                UTrueSize++;
        }
        size_t VminusUSize = n - UTrueSize;
        if (bipGraph.Y.size() >= (1.0 / 6) * VminusUSize) {
            for (auto v : bipGraph.X) {
                VC[v] = true;
                U[v] = true;
            }
            for (auto v : bipGraph.Y)
                U[v] = true;
        }

        if (bipGraph.Y.size() < (1.0 / 6) * VminusUSize || VminusUSize == 0)
            stop = true;
    }
}

Bipartite PlanarSeparatorVertexCover::bipartite(const NetworKit::Graph &G, std::vector<bool> &U,
                                                std::vector<int> &deg_residual) {
    std::queue<NetworKit::node> Q;
    std::vector<int> label(G.upperNodeIdBound(), 0);
    G.forNodes([&](NetworKit::node t) {
        if (!U[t] && deg_residual[t] > 2) {
            label[t] = '+';
            Q.push(t);
        }
    });
    while (!Q.empty()) {
        auto v = Q.front();
        Q.pop();

        G.forNeighborsOf(v, [&](NetworKit::node t) {
            if (!U[t] && label[t] == 0) {
                if (label[v] == '+')
                    label[t] = '-';
                else
                    label[t] = '+';

                Q.push(t);
            }
        });
    }

    std::vector<bool> isInX(G.upperNodeIdBound(), false);
    std::vector<bool> isInY(G.upperNodeIdBound(), false);
    std::vector<int> degB(G.upperNodeIdBound(), 0);

    for (int i = 0; i < static_cast<int>(label.size()); i++) {
        if (label[i] == '+') {
            isInX[i] = true;
            G.forNeighborsOf(i, [&](NetworKit::node t) {
                if (label[t] == '-') {
                    degB[i]++;
                    degB[t]++;
                }
            });
        }
        if (label[i] == '-')
            isInY[i] = true;
    }

    G.forNodes([&](NetworKit::node v) {
        if (!U[v] && degB[v] < 2)
            Q.push(v);
    });

    while (!Q.empty()) {
        auto v = Q.front();
        Q.pop();
        if (!isInX[v] && !isInY[v])
            continue;

        isInX[v] = false;
        isInY[v] = false;
        if (degB[v] == 0)
            continue;
        G.forNeighborsOf(v, [&](NetworKit::node t) {
            bool bEdge =
                (label[v] == '+' && label[t] == '-') || (label[v] == '-' && label[t] == '+');
            if (bEdge && (isInX[t] || isInY[t]) && --degB[t] < 2)
                Q.push(t);
        });
    }

    Bipartite bip = {{}, {}};

    for (size_t i = 0; i < isInX.size(); i++) {
        if (isInX[i])
            bip.X.push_back(i);
    }

    for (size_t i = 0; i < isInY.size(); i++) {
        if (isInY[i])
            bip.Y.push_back(i);
    }

    return bip;
}

std::vector<NetworKit::node> independentSetSolver(const NetworKit::Graph &cc) {
    std::vector<NetworKit::node> indpSet;
    std::vector<NetworKit::node> mapCompactToOriginal(cc.numberOfNodes());
    std::unordered_map<NetworKit::node, int> mapOriginalToCompact;
    int i = 0;
    cc.forNodes([&](NetworKit::node v) {
        mapOriginalToCompact[v] = i;
        mapCompactToOriginal[i++] = v;
    });
    std::vector<unsigned long long> adj(cc.numberOfNodes());
    cc.forEdges([&](NetworKit::node eu, NetworKit::node ev) {
        auto v = mapOriginalToCompact[ev];
        auto u = mapOriginalToCompact[eu];
        adj[v] = adj[v] | (1ULL << u);
        adj[u] = adj[u] | (1ULL << v);
    });

    unsigned long long maxIndpSet = 0;
    unsigned int maxSize = 0;
    for (size_t i = 0; i < (1ULL << adj.size()); i++) {
        bool valid = true;
        for (size_t j = 0; j < adj.size(); j++) {
            if (((i >> j) & 1ULL) && (adj[j] & i) != 0) {
                valid = false;
                break;
            }
        }
        if (valid) {
            size_t count = std::popcount(i);
            if (count > maxSize) {
                maxSize = count;
                maxIndpSet = i;
            }
        }
    }
    for (size_t i = 0; i < adj.size(); i++) {
        if ((maxIndpSet & (1ULL << i)) != 0)
            indpSet.push_back(mapCompactToOriginal[i]);
    }

    return indpSet;
}

void PlanarSeparatorVertexCover::run() {
    const auto &G = graph.value();
    std::vector<bool> VC(G.upperNodeIdBound());
    std::vector<bool> U(G.upperNodeIdBound());
    size_t n = G.upperNodeIdBound();
    prepare(G, U, VC, n);

    size_t UTrueSize = 0;
    for (auto v : U) {
        if (v)
            UTrueSize++;
    }
    if (n - UTrueSize == 0) {
        for (size_t i = 0; i < VC.size(); i++) {
            if (VC[i])
                vertex_cover.insert(i);
        }
    } else {
        std::unordered_set<NetworKit::node> residualNodes;
        G.forNodes([&](NetworKit::node v) {
            if (!U[v])
                residualNodes.insert(v);
        });
        const NetworKit::Graph residualGraph =
            NetworKit::GraphTools::subgraphFromNodes(G, residualNodes);
        double loglog = std::max(1.0, std::log2(std::log2((double)std::max<size_t>(4, n))));
        double epsilon = loglog / std::max<size_t>(1, residualGraph.numberOfNodes());
        MISP<NetworKit::node> mispAlgo(residualGraph, epsilon, independentSetSolver);
        mispAlgo.run();
        auto indepSet = mispAlgo.getIndependentSet();
        std::unordered_set<NetworKit::node> isSet(indepSet.begin(), indepSet.end());

        for (auto v : residualNodes)
            if (!isSet.count(v))
                vertex_cover.insert(v);
        for (size_t i = 0; i < VC.size(); i++)
            if (VC[i])
                vertex_cover.insert(i);
    }

    hasRun = true;
}
} // namespace Koala
