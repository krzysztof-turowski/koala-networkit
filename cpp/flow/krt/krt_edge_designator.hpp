#ifndef BETTER_MAX_FLOW_KRT_EDGE_DESIGNATOR_HPP
#define BETTER_MAX_FLOW_KRT_EDGE_DESIGNATOR_HPP

#include <bits/stdc++.h>

#include "edge_designator.hpp"

#ifndef PII
#define PII std::pair<int, int>
#endif
#ifndef LD
#define LD long double
#endif

std::vector<PII > get_undirected(std::vector<std::vector<KRTEdge>> const &g) {
    std::set<PII > edges;
    for (int i = 1; i < g.size(); ++i) {
        for (auto const &e: g[i]) {
            edges.insert(std::make_pair(i, e.dest));
            edges.insert(std::make_pair(e.dest, i));
        }
    }
    return std::vector<PII >(edges.begin(), edges.end());
}

class KRTEdgeDesignator : public EdgeDesignator {
    // Constraint: R0 * L/X >= 176
    // T = ceil(log(N) / log((R0 * L)/(88x))) + 4;
    // Parameters selected empirically to (1) match given constraints, (2) get best results on our datasets
    // TODO: Calculate parameters according to the Theorem for large enough graphs
    static constexpr LD R0 = 0.7, X = 2;
    static constexpr int L = 512;
    static constexpr int T = 7;

    int N, M, MAX_K;
    std::vector<int> degU, designated, rl, erl;
    std::vector<LD> ratios;
    std::unordered_set<int> U_prim, V_prim;
    std::vector<std::vector<int>> U, V;
    std::vector<std::vector<std::unordered_set<int>>> U_neighbors;

    void init_prim();

    void init_ratios();

    void init_neighbors();

    int deg(int) const;

    std::unordered_set<int> get_indexed_U(int);

    std::unordered_set<int> get_indexed_V(int);

    int encodeId(int, int) const;

    int decodeId(int) const;

    void update_rl(int v) {
        // Calculate r(v) - number of designated edges in U_prim x {v}
        int num_of_designated = 0;
//        for (auto u: U_prim) { // todo: it probably could be done faster
//            if (designated[u] == v) num_of_designated++;
//        }
        for (auto u: V[v]) {
            if (designated[u] == v && U_prim.count(u)) num_of_designated++;
        }
        LD r_v = static_cast<LD>(num_of_designated) / V[v].size();

        if (!V_prim.count(v) || r_v < R0) {
            rl[v] = 0;
        } else {
            for (int i = T - 1; i >= 0; --i) {
                if (r_v >= ratios[i]) {
                    rl[v] = i + 1;
                    return;
                }
            }
        }
    }

    void update_erl(int v) {
        erl[v] = rl[v];

        for (auto u: V[v]) {
            if (!U_prim.count(u)) continue;
            bool found = false;
            for (int ratio = 0; ratio <= T; ++ratio) {
                if (U_neighbors[u][ratio].count(v)) {
                    found = true;
                    U_neighbors[u][ratio].erase(v);
                    break;
                }
            }
            if (found) U_neighbors[u][erl[v]].insert(v);
        }
    }

    void remove_edge(int u, int v) {
        // Remove v from u's neighbor lists
        for (int ratio = 0; ratio <= T; ++ratio) {
            if (U_neighbors[u][ratio].count(v)) {
                U_neighbors[u][ratio].erase(v);
                degU[u]--;
                break;
            }
        }

//        org version
//        if (U_prim.count(u) && deg(u) < L) {
//            U_prim.erase(u);
//            if (designated[u] == v) {
//                update_rl(v);
//                if (rl[v] < erl[v] - 1) {
//                    update_erl(v);
//                }
//            }
//        }

        if (designated[u] == v) {
            designated[u] = -1;
            if (U_prim.count(u) && deg(u) < L) {
                U_prim.erase(u);
            }
            if (U_prim.count(u)) {
                update_rl(v);
                if (rl[v] < erl[v] - 1) {
                    update_erl(v);
                }
            }
            designate_edge(u);
        }
    }

    int designate_edge(int u) {
        int v = -1;
        if (!U_prim.count(u)) {
            // u in U \ U' => designate any incident edge
            // Czy tu nie powinien byc lookup do pierwszej niepustej listy? nie musi bo dla takich u jak w zalozzeniu U_nei sie nigdy nie zmienia
            if (!U_neighbors[u].front().empty()) v = *(U_neighbors[u].front().begin());
            designated[u] = v;
        } else {
            for (int ratio = 0; ratio <= T; ++ratio) {
                if (!U_neighbors[u][ratio].empty()) {
                    v = *(U_neighbors[u][ratio].begin());
                    break;
                }
            }

            designated[u] = v;

            update_rl(v);

            if (rl[v] > erl[v]) {
                update_erl(v);
            }
        }
        return v;
    }

    LD reset() {
        auto k = T;

        while (get_indexed_U(k - 3).size() >= (ratios[k - 3] * L) * get_indexed_U(k).size() / (88.0 * X)) {
            k -= 3;
        }

        auto WK1 = get_indexed_U(k - 1);
        auto VK1 = get_indexed_V(k - 1);
        for (auto v: VK1) {
            while (rl[v] >= k - 1) {
                for (auto u: U_prim) {
                    if (designated[u] == v) {
                        designated[u] = -1;
                        break; // ?
                    }
                }
                update_rl(v);
            }
            if (rl[v] > erl[v] || rl[v] < erl[v] - 1) {
                update_erl(v);
            }
        }
        for (auto u: WK1) {
            if (designated[u] == -1) {
                designate_edge(u);
            }
        }

        return k;
    }


public:

    KRTEdgeDesignator(int n, const std::vector<std::vector<KRTEdge>> &g);

    int ce(int, int) override;

    void response_adversary(int a, int da, int b, int db, bool remove_v) override {
        if (!remove_v) {
            int u = encodeId(a, da);
            int v = encodeId(b, db);

            remove_edge(u, v);
            if (designated[u] == v) {
                auto v_prim = designate_edge(u);
                if (rl[v_prim] == T) {
                    bool all_less;
                    do {
                        reset();
                        all_less = true;

                        for (int i = 0; i < N; ++i) {
                            if (rl[i] >= T) {
                                all_less = false;
                                break;
                            }
                        }
                    } while (all_less);
                }
            }
        } else {
            int v = encodeId(b, db);
            // respond as if the adversary removed each edge (u, v) in any sequence
            for (auto nei: V[v]) {
                remove_edge(nei, v);
            }
        }
    }
};

void KRTEdgeDesignator::init_prim() {
    U_prim.reserve(2 * N);
    V_prim.reserve(2 * N);
    for (int i = 0; i < N; ++i) {
        if (U[i].size() >= L) {
            U_prim.insert(i);
        }
        if (V[i].size() >= L) {
            V_prim.insert(i);
        }
    }
}

void KRTEdgeDesignator::init_ratios() {
    ratios = std::vector<LD>(T);
    ratios[0] = R0;
    for (int i = 1; i < T; ++i) {
        ratios[i] = (1 + static_cast<long double>(1.0) / X) * ratios[i - 1];
//        std::cout << ratios[i] << std::endl;
    }
}

void KRTEdgeDesignator::init_neighbors() {
    U_neighbors.clear();
    U_neighbors.resize(N);

    for (int i = 0; i < N; ++i) {
        U_neighbors[i].resize(T + 1);
        U_neighbors[i][0] = std::unordered_set<int>(U[i].begin(), U[i].end());
        degU[i] = static_cast<int>(U[i].size());
    }
}

int KRTEdgeDesignator::deg(int u) const {
    return degU[u];
    int deg = 0;
    for (int ratio = 0; ratio <= T; ++ratio) {
        deg += static_cast<int>(U_neighbors[u][ratio].size());
    }
    return deg;
}

std::unordered_set<int> KRTEdgeDesignator::get_indexed_U(int k) {
    std::unordered_set<int> UK;
    for (auto const u: U_prim) {
        int des = designated[u];
        if (V_prim.count(des) && rl[des] >= k) {
            UK.insert(u);
        }
    }
    return UK;
}

std::unordered_set<int> KRTEdgeDesignator::get_indexed_V(int k) {
    std::unordered_set<int> VK;
    for (auto const v: V_prim) {
        if (rl[v] >= k) {
            VK.insert(v);
        }
    }
    return VK;
}

KRTEdgeDesignator::KRTEdgeDesignator(int n, const std::vector<std::vector<KRTEdge>> &g) {
    auto A_undirected = get_undirected(g);

    MAX_K = 2 * n;
    N = (n + 1) * MAX_K;
    M = 0;

    U = std::vector<std::vector<int>>(N);
    V = std::vector<std::vector<int>>(N);
    degU = std::vector<int>(N, 0);
    designated = std::vector<int>(N, -1);
    rl = std::vector<int>(N, 0);
    erl = std::vector<int>(N, 0);

    for (auto const &uv: A_undirected) {
        int u = uv.first, v = uv.second;
        for (int k = 1; k < MAX_K; ++k) {
            int left = encodeId(u, k);
            int right = encodeId(v, k - 1);
            U[left].push_back(right);
            V[right].push_back(left);
            M++;
        }
    }
    init_ratios();
    init_neighbors();
    init_prim();

    for (int i = 0; i < N; ++i) {
        designate_edge(i);
    }
}

int KRTEdgeDesignator::ce(int i, int k) {
    return decodeId(designated[encodeId(i, k)]);
}

int KRTEdgeDesignator::encodeId(int i, int k) const {
    return i * MAX_K + k;
}

int KRTEdgeDesignator::decodeId(int i) const {
    if (i == -1) return -1;
    return i / MAX_K;
}

#endif //BETTER_MAX_FLOW_KRT_EDGE_DESIGNATOR_HPP
