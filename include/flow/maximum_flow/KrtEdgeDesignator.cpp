/*
 * KrtEdgeDesignator.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <map>
#include <vector>

namespace Koala {

void KRTEdgeDesignator::initialize_prim() {
    U_prim.reserve(2 * N), V_prim.reserve(2 * N);
    for (int i = 0; i < N; ++i) {
        if (U[i].size() >= L) {
            U_prim.insert(i);
        }
        if (V[i].size() >= L) {
            V_prim.insert(i);
        }
    }
}

void KRTEdgeDesignator::initialize_ratios() {
    ratios = std::vector<long double>(T);
    ratios[0] = R0;
    for (int i = 1; i < T; ++i) {
        ratios[i] = (1 + static_cast<long double>(1.0) / X) * ratios[i - 1];
    }
}

void KRTEdgeDesignator::initialize_neighbors() {
    U_neighbors.clear();
    U_neighbors.resize(N);

    for (int i = 0; i < N; ++i) {
        U_neighbors[i].resize(T + 1);
        U_neighbors[i][0] = std::unordered_set<int>(U[i].begin(), U[i].end());
        degU[i] = static_cast<int>(U[i].size());
    }
} 

int KRTEdgeDesignator::deg(int u) const{
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

int KRTEdgeDesignator::encodeId(int i, int k) const {
    return i * MAX_K + k;
}

int KRTEdgeDesignator::decodeId(int i) const {
    if (i == -1) return -1;
    return i / MAX_K;
}

void KRTEdgeDesignator::update_rl(int v) {
    // Calculate r(v) - number of designated edges in U_prim x {v}
    int num_of_designated = 0;
//        for (auto u: U_prim) { // todo: it probably could be done faster
//            if (designated[u] == v) num_of_designated++;
//        }
    for (auto u: V[v]) {
        if (designated[u] == v && U_prim.count(u)) num_of_designated++;
    }
    long double r_v = static_cast<long double>(num_of_designated) / V[v].size();

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

void KRTEdgeDesignator::update_erl(int v) {
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

void KRTEdgeDesignator::remove_edge(int u, int v) {
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

int KRTEdgeDesignator::designate_edge(int u) {
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

long double KRTEdgeDesignator::reset() {
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

void KRTEdgeDesignator::initialize(const std::optional<NetworKit::Graph> &graph) {
    int n = graph->numberOfNodes();
    MAX_K = 2 * n, N = (n + 1) * MAX_K, M = 0;
    U = std::vector<std::vector<int>>(N), V = std::vector<std::vector<int>>(N);
    degU = std::vector<int>(N, 0);
    designated = std::vector<int>(N, -1);
    rl = std::vector<int>(N, 0), erl = std::vector<int>(N, 0);

    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        for (int k = 1; k < MAX_K; ++k) {
            int left = encodeId(u + 1, k);
            int right = encodeId(v + 1, k - 1);
            U[left].push_back(right);
            V[right].push_back(left);
            M++;
        }
    });
    initialize_ratios();
    initialize_neighbors();
    initialize_prim();
    for (int i = 0; i < N; ++i) {
        designate_edge(i);
    }
}

int KRTEdgeDesignator::current_edge(int i, int k) {
    return decodeId(designated[encodeId(i, k)]);
}

void KRTEdgeDesignator::response_adversary(int a, int da, int b, int db, bool remove_v) {
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

} /* namespace Koala */
