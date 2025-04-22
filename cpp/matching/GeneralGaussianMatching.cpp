#include <matching/GeneralGaussianMatching.hpp>
#include <matching/GaussElimination.hpp>
#include <matching/utils.hpp>
#include <matching/BipartiteGaussianMatching.hpp>

#include <set>
#include <map>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace NetworKit;

namespace Koala {
    struct Partition {
        vector<Graph> subgraphs;
        vector<pair<int, int>> componentMap; // v -> c,i
        vector<vector<int>> components;

        int getLabel(int c, int i) {
            return components[c][i];
        }

        pair<int, int> getComponent(int v) {
            return componentMap[v];
        }
    };

    Partition partitionGraph(const Graph& graph, vector<set<int>>& components);
    vector<set<int>> getSimpleComponents(const Graph& graph, const MatrixXd& AG);
    vector<set<int>> getConnectedComponents(const Graph& graph, vector<bool> visited);
    set<int> getNontrivialClass(const MatrixXd& AG);

    Matching greedyMatching(const Graph& G);
    Matching generalMatching(const Graph& G, const MatrixXd& AG);
    Matching partition(const Graph& graph, const MatrixXd& AG);
    Matching simplePartition(const Graph& graph, const MatrixXd& AG);
    MatrixXd generateMatrix(const Graph& G);

    GeneralGaussianMatching::GeneralGaussianMatching(const Graph& G) : G(G) {};

    void GeneralGaussianMatching::run() {
        MatrixXd AG = generateMatrix(G);
        M = generalMatching(G, AG);
    }

    Matching GeneralGaussianMatching::getMatching() {
        return M;
    }


    Matching greedyMatching(const Graph& G) {
        int n = G.numberOfNodes();
        Matching M;

        vector<bool> isMatched(n, false);
        G.forEdges([&](node u, node v) {
            if (!isMatched[u] && !isMatched[v]) {
                isMatched[u] = isMatched[v] = true;
                M.insert({ u,v });
            }
            });

        return M;
    }

    Matching getMaximalMatching(const MatrixXd& AG, const Matching& M) {
        int n = AG.cols();
        int k = 2 * M.size();

        vector<int> rowIndices(n), colIndices(n);

        int i = 0;
        for (auto [u, v] : M) {
            rowIndices[i] = colIndices[i + 1] = u;
            rowIndices[i + 1] = colIndices[i] = v;
            i += 2;
        }
        for (; i < n; ++i) {
            rowIndices[i] = colIndices[i] = i;
        }

        MatrixXd AGr(n, n);
        for (int i = 0; i < n; ++i) {
            AGr.row(i) = AG.row(rowIndices[i]);
        }

        MatrixXd AGrc(n, n);
        for (int i = 0; i < n; ++i) {
            AGrc.col(i) = AGr.col(colIndices[i]);
        }

        auto eliminated = simpleElimination(AGrc, k);
        Matching M1;
        for (auto i : eliminated) {
            int u = min(rowIndices[i], colIndices[i]);
            int v = max(rowIndices[i], colIndices[i]);
            M1.insert({ u,v });
        }
        return M1;
    }

    Matching partition(const Graph& G, const MatrixXd& AG) {
        int n = AG.cols();

        auto components = getSimpleComponents(G, AG);
        auto partition = partitionGraph(G, components);

        Matching M;
        for (int c = 0; c < partition.subgraphs.size(); ++c) {
            auto Gc = partition.subgraphs[c];
            auto AGc = generateMatrix(Gc);
            auto M1 = simplePartition(Gc, AGc);
            for (auto [ui, vi] : M1) {
                int u = partition.getLabel(c, ui);
                int v = partition.getLabel(c, vi);
                M.insert({ u,v });
            }
        }

        return M;
    }

    Matching simplePartition(const Graph& G, const MatrixXd& AG) {
        auto Sv = getNontrivialClass(AG);
        if (Sv.size() == 0) {
            return greedyMatching(G); // TODO
        }
        vector<bool> isInS(G.numberOfNodes(), false);
        for (auto v : Sv) {
            isInS[v] = true;
        }

        // Components of C1...Ck of G-S and S
        auto components = getConnectedComponents(G, isInS);
        components.push_back(Sv);
        auto partition = partitionGraph(G, components);

        auto S = partition.subgraphs.back();
        int sSize = S.numberOfNodes();
        assert(sSize == partition.subgraphs.size() - 1); // |{C1, C2, ...}| = sSize = |S|

        // Contract Ci to single vertex
        S.addNodes(sSize);
        G.forEdges([&](node u, node v) {
            if (isInS[u] && isInS[v]) {
                auto [cu, ui] = partition.getComponent(u);
                auto [cv, vi] = partition.getComponent(v);
                S.removeEdge(ui, vi);
            } else if (isInS[u] != isInS[v]) {
                int w = isInS[u] ? v : u;
                int s = isInS[u] ? u : v;

                auto [cw, wi] = partition.getComponent(w);
                auto [sw, si] = partition.getComponent(s);
                S.addEdge(si, sSize + cw);
            }
            });

        // Get the matching of S and contracted ci
        auto bpMatch = BipartiteGaussianMatching(S);
        bpMatch.run();
        auto matchS = bpMatch.getMatching();

        Matching M;
        for (auto [si, ci] : matchS) {
            // si has lower index than a contracted vertex ci
            if (si > ci) swap(si, ci);
            int s = partition.getLabel(sSize, si);
            int c = ci - sSize;

            // get some neighbor of s from component c
            int v = -1;
            G.forEdgesOf(s, [&](node u) {
                if (partition.getComponent(u).first == c) {
                    v = u;
                }
                });
            assert(v != -1);

            M.insert({ s,v });

            auto Gc = partition.subgraphs[c];
            if (Gc.numberOfNodes() == 1)
                continue;

            // Swap v and the last node for quick delete
            int vi = partition.getComponent(v).second;
            int ui = partition.components[c].size() - 1;
            int u = partition.getLabel(c, ui);

            Gc.forEdgesOf(vi, [&](int wi) {
                Gc.removeEdge(vi, wi);
                });
            Gc.forEdgesOf(ui, [&](int wi) {
                if (vi == wi) return;
                Gc.addEdge(vi, wi);
                });
            Gc.removeNode(ui);

            partition.components[c][vi] = u;
            partition.components[c].pop_back();
            partition.componentMap[u] = { c, vi };
            partition.componentMap[v] = { -1, -1 };

            // Match the rest of component c
            auto genMatch = GeneralGaussianMatching(Gc);
            genMatch.run();
            auto matchC = genMatch.getMatching();

            for (auto [ui, vi] : matchC) {
                int u = partition.getLabel(c, ui);
                int v = partition.getLabel(c, vi);
                M.insert({ u,v });
            }
        }

        return M;
    }

    Matching generalMatching(const Graph& G, const MatrixXd& AG) {
        int n = G.numberOfNodes();
        if (n == 0) return {};

        Matching M = greedyMatching(G);
        Matching M1 = getMaximalMatching(AG, M);

        vector<set<int>> components(2, set<int>());
        for (auto [u, v] : M1) {
            components[0].insert(u);
            components[0].insert(v);
        }
        for (int v = 0; v < n; ++v) {
            if (components[0].find(v) == components[0].end()) {
                components[1].insert(v);
            }
        }

        auto part = partitionGraph(G, components);
        auto G1 = part.subgraphs[1];
        auto AG1 = generateMatrix(G1);

        Matching M2;
        if (8 * M1.size() >= n) {
            M2 = generalMatching(G1, AG1);
        } else {
            M2 = partition(G1, AG1);
        }

        for (auto [ui, vi] : M2) {
            M1.insert({ part.getLabel(1, ui), part.getLabel(1, vi) });
        }

        return M1;
    }


    // TODO: move to some utils
    void dfs(int u, vector<bool>& visited, set<int>& connected, function<vector<int>(int)>& edgesOf) {
        visited[u] = true;
        connected.insert(u);
        for (auto v : edgesOf(u)) {
            if (visited[v]) continue;
            dfs(v, visited, connected, edgesOf);
        }
    };

    vector<set<int>> getSimpleComponents(const Graph& G, const MatrixXd& AG) {
        int n = G.numberOfNodes();

        vector<bool> visited(n, false);
        vector<set<int>> connected;
        function<vector<int>(int)> edgesOf = [&](int u) {
            vector<int> edges;
            for (int v = 0; v < n; ++v) {
                if (G.hasEdge(u, v) && !eq(AG(u, v), 0)) {
                    edges.push_back(v);
                }
            }
            return edges;
            };

        for (int v = 0; v < n; ++v) {
            if (visited[v]) continue;
            connected.push_back(set<int>());
            dfs(v, visited, connected[connected.size() - 1], edgesOf);
        }
        return connected;
    }

    vector<set<int>> getConnectedComponents(const Graph& G, vector<bool> visited) {
        int n = G.numberOfNodes();

        vector<set<int>> connected;
        function<vector<int>(int)> edgesOf = [&](int u) {
            vector<int> edges;
            G.forEdgesOf(u, [&](int v) {
                edges.push_back(v);
                });
            return edges;
            };

        for (int v = 0; v < n; ++v) {
            if (visited[v]) continue;
            connected.push_back(set<int>());
            dfs(v, visited, connected[connected.size() - 1], edgesOf);
        }
        return connected;
    }

    set<int> getNontrivialClass(const MatrixXd& AG) {
        int n = AG.cols();

        vector<bool> visited(n, false);
        set<int> S;
        function<vector<int>(int)> edgesOf = [&](int u) {
            vector<int> edges;
            for (int v = 0; v < n; ++v) {
                if (eq(AG(u, v), 0)) {
                    edges.push_back(v);
                }
            }
            return edges;
            };

        for (int u = 0; u < n; ++u) {
            for (int v = u + 1; v < n; ++v) {
                if (eq(AG(u, v), 0)) {
                    dfs(u, visited, S, edgesOf);
                    return S;
                }
            }
        }
        return {};
    }

    Partition partitionGraph(const Graph& G, vector<set<int>>& components) {
        int n = G.numberOfNodes();
        int m = components.size();

        Partition partition;
        partition.subgraphs.resize(m);
        partition.componentMap.resize(n);
        partition.components.resize(m);

        for (int c = 0; c < m; ++c) {
            int nc = components[c].size();
            partition.components[c].resize(nc);

            int i = 0;
            for (auto& v : components[c]) {
                partition.componentMap[v] = { c, i };
                partition.components[c][i] = v;
                i++;
            }

            partition.subgraphs[c] = Graph(nc, false, false);
        }

        // copy G
        G.forEdges([&](node u, node v) {
            auto [cu, ui] = partition.componentMap[u];
            auto [cv, vi] = partition.componentMap[v];
            if (cu == cv) {
                partition.subgraphs[cu].addEdge(ui, vi);
            }
            });

        return partition;
    }

    MatrixXd generateMatrix(const Graph& G) {
        int n = G.numberOfNodes();

        MatrixXd AG = ArrayXXd::Zero(n, n);
        G.forEdges([&](node u, node v) {
            double Xuv = generateRandom();
            AG(u, v) = Xuv;
            AG(v, u) = -Xuv;
            });
        return AG.inverse();
    }
}

