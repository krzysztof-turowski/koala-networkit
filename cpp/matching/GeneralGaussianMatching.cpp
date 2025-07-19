#include <matching/GeneralGaussianMatching.hpp>
#include <matching/GaussElimination.hpp>
#include <matching/utils.hpp>
#include <matching/BipartiteGaussianMatching.hpp>
#include <matching/DynamicComponents.hpp>

#include <set>
#include <map>

#include <Eigen/Dense>
#include <networkit/graph/GraphTools.hpp>

using namespace std;
using namespace Eigen;
using namespace NetworKit;

namespace Koala {
    void dfs(int u, vector<bool>& visited, set<int>& connected, function<vector<int>(int)> edgesOf);
    vector<set<int>> getSimpleComponents(const Graph& graph, const MatrixXd& AG);
    vector<set<int>> getConnectedComponents(const Graph& graph, vector<bool> visited);
    set<int> getNontrivialClass(const MatrixXd& AG);

    Matching greedyMatching(const Graph& G);
    Matching generalMatching(const Graph& G, const MatrixXd& AG);
    Matching partition(const Graph& graph, const MatrixXd& AG);
    Matching simplePartition(const Graph& graph, const MatrixXd& AG);
    MatrixXd generateMatrix(const Graph& G);

    GeneralGaussianMatching::GeneralGaussianMatching(const Graph& G1) {
        if (G1.numberOfNodes() == 0) return;
        auto [_G, _oldIdx] = reindexGraph(G1);
        G = _G;
        oldIdx = _oldIdx;
    };

    void GeneralGaussianMatching::run() {
        if (G.numberOfNodes() == 0) return;
        MatrixXd AG = generateMatrix(G);
        M = generalMatching(G, AG);
    }

    Matching GeneralGaussianMatching::getMatching() {
        Matching M1;
        for (auto [u, v] : M) {
            M1.insert({ oldIdx[u], oldIdx[v] });
        }
        return M1;
    }


    Matching greedyMatching(const Graph& G) {
        int n = G.numberOfNodes();
        Matching M;

        vector<bool> isMatched(n, false);
        for (const auto& [u, v] : G.edgeRange()) {
            if (!isMatched[u] && !isMatched[v]) {
                isMatched[u] = isMatched[v] = true;
                M.insert({ u,v });
            }
        }
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
        cout << "PARTITION G:" << G.numberOfNodes() << endl;
        int n = G.numberOfNodes();

        auto components = getSimpleComponents(G, AG);
        vector<tuple<Graph, vector<int>, MatrixXd>> subgraphs(components.size());
        for (int i = 0; i < components.size(); ++i) {
            subgraphs[i] = reindexGraph(GraphTools::subgraphFromNodes(G, components[i].begin(), components[i].end()), AG);
        }

        Matching M;
        cout << "Partition of " << subgraphs.size() << " simple components" << endl;
        for (int c = 0; c < subgraphs.size(); ++c) {
            auto [Gc, oldIdx, AGc] = subgraphs[c];
            auto M1 = simplePartition(Gc, AGc);
            // auto AGc = generateMatrix(Gc);

            assert(M1.size() == Gc.numberOfNodes() / 2);
            for (auto [u1, v1] : M1) {
                int u = oldIdx[u1], v = oldIdx[v1];
                M.insert({ u, v });

                assert(u != v);
                assert(components[c].find(u) != components[c].end());
                assert(components[c].find(v) != components[c].end());
            }
        }
        assert(M.size() == G.numberOfNodes() / 2);
        return M;
    }

    Matching simplePartition(const Graph& G, const MatrixXd& AG) {
        cout << "SIMPLE_PARTITION G:" << endl;

        Matching M;
        DynamicComponents DC(G);

        auto Sv = getNontrivialClass(AG);
        if (Sv.size() == 0) {
            GeneralGaussianMatching gen(G);
            gen.run();
            M = gen.getMatching();
            assert(M.size() == G.numberOfNodes() / 2);
            return M;
        }
        int sSize = Sv.size();
        vector<bool> isInS(G.numberOfNodes(), false);
        for (auto v : Sv) {
            isInS[v] = true;
        }

        set<int> Tv;
        for (auto sv : Sv) {
            for (const auto& t : G.neighborRange(sv)) {
                DC.removeEdge(sv, t);
                if (!isInS[t]) Tv.insert(t);
            }
        }

        // Get largest component of G\S
        int maxComponentSize = 0, maxComponentV;
        for (auto t : Tv) {
            int s = DC.getComponentSize(t);
            if (s > maxComponentSize) {
                maxComponentSize = s;
                maxComponentV = t;
            }
        }

        // Get vertices of smaller component that were connected to S
        set<int> Dv;
        for (auto u : Tv) {
            if (!DC.isConnected(maxComponentV, u)) {
                Dv.insert(u);
            }
        }

        // Get vertices of all components of G\S and S
        vector<set<int>> components({ {} });
        for (auto v : G.nodeRange()) components[0].insert(v);
        for (auto v : Sv) components[0].erase(v);

        for (auto v : Dv) {
            components.push_back(set<int>({}));
            auto& comp = components[components.size() - 1];

            vector<bool> visited(G.numberOfNodes(), false);
            dfs(v, visited, comp, [&](int u) {
                auto it = DC.G.neighborRange(u);
                return vector<int>(it.begin(), it.end());
                });
            for (auto v : comp) components[0].erase(v);
        }

        int s = DC.G.addNode();
        for (auto v : Tv) {
            DC.addEdge(s, v);
        }
        for (auto& comp : components) {
            comp.insert(s);
        }
        components.push_back(Sv);

        vector<int> componentOf(G.numberOfNodes());
        for (int i = 0; auto& comp :components) {
            for (auto v : comp) {
                if (v == G.numberOfNodes()) continue;
                componentOf[v] = i;
            }
            i++;
        }

        vector<Graph> subgraphs(components.size());
        for (int i = 0; i < components.size(); ++i) {
            subgraphs[i] = GraphTools::subgraphFromNodes(DC.G, components[i].begin(), components[i].end());
        }
        auto& SG = subgraphs[sSize];
        auto& C1G = subgraphs[0];

        // Match the biggest component
        cout << "Match C0 component" << endl;
        GeneralGaussianMatching genC1(C1G); // TODO non-trivial class
        genC1.run();
        auto MC1 = genC1.getMatching();
        assert(MC1.size() == C1G.numberOfNodes() / 2);
        for (auto [u, v] : MC1) assert(u != v);


        // Remove vertex matched with the contracted biggest component
        cout << "Removed s matched with C0" << endl;
        for (auto [u, v] : MC1) {
            cout << u << ' ' << v << ' ' << s << endl;
            if (u == s || v == s) {
                int cv = v == s ? u : v;
                for (auto sv : Sv) {
                    if (G.hasEdge(cv, sv)) {
                        assert(sv != G.numberOfNodes());
                        assert(cv != G.numberOfNodes());
                        M.insert({ cv, sv });
                        SG.removeNode(sv);
                        Sv.erase(sv);
                        break;
                    }
                }
            } else {
                assert(u != G.numberOfNodes());
                assert(v != G.numberOfNodes());
                M.insert({ u,v });
            }
        }
        assert(Sv.size() == sSize - 1);

        // Create a graph of S and a contracted vertex for each of the smaller components of G\S
        cout << "contract smaller components" << endl;
        vector<int> contractedNodes(sSize);
        map<int, int> contractedComponents;
        for (int i = 1; i < sSize; ++i) {
            int v = SG.addNode();
            contractedNodes[i] = v;
            contractedComponents[v] = i;
        }
        for (auto v : Dv) {
            for (int i = 0; auto sv : Sv) {
                int cv = componentOf[v];
                if (cv == 0) continue;
                if (G.hasEdge(sv, v)) {
                    SG.addEdge(sv, contractedNodes[cv]);
                }
            }
        }

        cout << "bp matching of S:" << endl;
        BipartiteGaussianMatching bpS(SG);
        bpS.run();
        auto MS = bpS.getMatching();
        assert(MS.size() == SG.numberOfNodes() / 2);
        for (auto [u, v] : MS) assert(u != v);

        // Eliminated matched vertex s with some vertex from corresponding smaller component
        cout << "Add matching from S:" << endl;
        for (auto [s, cv] : MS) {
            if (contractedComponents.find(cv) == contractedComponents.end())
                swap(s, cv);

            int c = contractedComponents[cv];
        }
        for (auto [s, cv] : MS) {
            if (contractedComponents.find(cv) == contractedComponents.end())
                swap(s, cv);
            assert(contractedComponents.find(cv) != contractedComponents.end());

            int c = contractedComponents[cv];
            for (auto v : components[c]) {
                if (G.hasEdge(s, v)) {
                    M.insert({ s, v });
                    subgraphs[c].removeNode(v);
                    subgraphs[c].removeNode(s);
                    components[c].erase(v);
                    break;
                }
                assert(false);
            }
        }

        // Match smaller components
        for (int i = 1; i < sSize; ++i) {
            auto& SCi = subgraphs[i];
            cout << "Match C" << i << endl;
            GeneralGaussianMatching gen(SCi);
            gen.run();
            auto MCi = gen.getMatching();
            assert(MCi.size() == SCi.numberOfNodes() / 2);
            for (auto [u, v] : MCi) {
                assert(u != v);
                assert(u != G.numberOfNodes());
                assert(v != G.numberOfNodes());
                M.insert({ u,v });
            }
        }

        assert(M.size() == G.numberOfNodes() / 2);
        return M;
    }

    Matching generalMatching(const Graph& G, const MatrixXd& AG) {
        cout << "GENERAL_MATCHING G:" << G.numberOfNodes() << endl;

        const int n = G.numberOfNodes();
        if (n == 0) return {};

        Matching M = greedyMatching(G);
        if (2 * M.size() == n) {
            return M;
        }
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

        auto [G1, oldIdx] = reindexGraph(GraphTools::subgraphFromNodes(G, components[1].begin(), components[1].end()));
        auto AG1 = generateMatrix(G1);

        Matching M2;
        if (8 * M1.size() >= n) {
            M2 = generalMatching(G1, AG1);
        } else {
            M2 = partition(G1, AG1);
        }
        assert(M2.size() == G1.numberOfNodes() / 2);
        for (auto [u, v] : M2) {
            assert(u != v);
            M1.insert({ oldIdx[u], oldIdx[v] });
        }

        return M1;
    }

    void dfs(int u, vector<bool>& visited, set<int>& connected, function<vector<int>(int)> edgesOf) {
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
            for (const auto& v : G.neighborRange(u)) {
                if (!eq(AG(u, v), 0)) {
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
        auto edgesOf = [&](int u) {
            return vector<int>(G.neighborRange(u).begin(), G.neighborRange(u).end());
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

    MatrixXd generateMatrix(const Graph& G) {
        int n = G.numberOfNodes();

        MatrixXd AG = ArrayXXd::Zero(n, n);
        G.forEdges([&](node u, node v) {
            auto Xuv = generateRandom();
            AG(u, v) = Xuv;
            AG(v, u) = -Xuv;
            });
        return AG.inverse();
    }
}

