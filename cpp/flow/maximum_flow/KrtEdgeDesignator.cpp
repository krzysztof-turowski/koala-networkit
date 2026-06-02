/*
 * KrtEdgeDesignator.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "flow/maximum_flow/KrtEdgeDesignator.hpp"

namespace Koala {

void KRTEdgeDesignator::initialize_prim() {
    U_prim.clear(), V_prim.clear();
    U_prim.reserve(2 * N), V_prim.reserve(2 * N);
    for (NetworKit::node i = 0; i < N; ++i) {
        if (U[i].size() >= l) {
            U_prim.insert(i);
        }
        if (V[i].size() >= l) {
            V_prim.insert(i);
        }
    }
}

void KRTEdgeDesignator::initialize_parameters(const Parameters &parameters) {
    r0 = parameters.R0.value_or(0.7);
    x = parameters.X.value_or(2);
    if (r0 <= 0 || x < 2) {
        throw std::invalid_argument("KRTEdgeDesignator requires r0 > 0 and x >= 2");
    }

    l = parameters.L.value_or(
        static_cast<NetworKit::count>(std::floor(176 * x / r0)) + 1);
    if (l == 0 || r0 * l / x <= 176) {
        throw std::invalid_argument("KRTEdgeDesignator requires r0 * l / x > 176");
    }

    const long double denominator = std::log(r0 * l / (88 * x));
    t = parameters.T.value_or(
        3 * static_cast<int>(std::ceil(std::log(std::max<NetworKit::count>(N, 2))
                                      / denominator))
            + 4);
    if (t < 4) {
        throw std::invalid_argument("KRTEdgeDesignator requires t >= 4");
    }
}

void KRTEdgeDesignator::initialize_ratios() {
    ratios = std::vector<long double>(t);
    ratios[0] = r0;
    for (int i = 1; i < t; ++i) {
        ratios[i] = (1 + static_cast<long double>(1.0) / x) * ratios[i - 1];
    }
}

void KRTEdgeDesignator::initialize_neighbors() {
    U_neighbors.clear();
    U_neighbors.resize(N);
    for (NetworKit::node i = 0; i < N; ++i) {
        U_neighbors[i].resize(t + 1);
        U_neighbors[i][0] =
            std::unordered_set<NetworKit::node>(U[i].begin(), U[i].end());
        deg_U[i] = U[i].size();
    }
}

std::unordered_set<NetworKit::node> KRTEdgeDesignator::get_indexed_U(int k) {
    std::unordered_set<NetworKit::node> UK;
    for (auto const u : U_prim) {
        NetworKit::node des = designated[u];
        if (V_prim.count(des) && rl[des] >= k) {
            UK.insert(u);
        }
    }
    return UK;
}

std::unordered_set<NetworKit::node> KRTEdgeDesignator::get_indexed_V(int k) {
    std::unordered_set<NetworKit::node> VK;
    for (auto const v : V_prim) {
        if (rl[v] >= k) {
            VK.insert(v);
        }
    }
    return VK;
}

NetworKit::node KRTEdgeDesignator::encode_id(NetworKit::node i, int k) const {
    return i * MAX_K + k;
}

NetworKit::node KRTEdgeDesignator::decode_id(NetworKit::node i) const {
    return i != NetworKit::none ? i / MAX_K : NetworKit::none;
}

void KRTEdgeDesignator::update_rl(NetworKit::node v) {
    if (v == NetworKit::none || V[v].empty()) {
        return;
    }
    NetworKit::count num_of_designated = 0;
    for (auto u : V[v]) {
        if (designated[u] == v && U_prim.count(u)) {
            num_of_designated++;
        }
    }
    long double r_v = static_cast<long double>(num_of_designated) / V[v].size();

    if (!V_prim.count(v) || r_v < r0) {
        rl[v] = 0;
    } else {
        for (int i = t - 1; i >= 0; --i) {
            if (r_v >= ratios[i]) {
                rl[v] = i + 1;
                return;
            }
        }
    }
}

void KRTEdgeDesignator::update_erl(NetworKit::node v) {
    erl[v] = rl[v];
    for (auto u : V[v]) {
        if (!U_prim.count(u)) continue;
        bool found = false;
        for (int ratio = 0; ratio <= t; ++ratio) {
            if (U_neighbors[u][ratio].count(v)) {
                found = true;
                U_neighbors[u][ratio].erase(v);
                break;
            }
        }
        if (found) {
            U_neighbors[u][erl[v]].insert(v);
        }
    }
}

bool KRTEdgeDesignator::remove_edge(NetworKit::node u, NetworKit::node v) {
    bool found = false;
    for (int ratio = 0; ratio <= t; ++ratio) {
        if (U_neighbors[u][ratio].count(v)) {
            U_neighbors[u][ratio].erase(v);
            deg_U[u]--;
            found = true;
            break;
        }
    }
    if (!found) {
        return false;
    }

    const bool was_designated = designated[u] == v;
    const bool was_in_U_prim = U_prim.count(u);
    const NetworKit::node old_designation = designated[u];
    if (was_designated) {
        designated[u] = NetworKit::none;
    }

    const bool shifted = was_in_U_prim && deg_U[u] < l;
    if (shifted) {
        U_prim.erase(u);
    }
    if (was_in_U_prim && (shifted || was_designated)) {
        update_rl(old_designation);
        if (old_designation != NetworKit::none
                && rl[old_designation] < erl[old_designation] - 1) {
            update_erl(old_designation);
        }
    }

    if (shifted) {
        for (int ratio = 1; ratio <= t; ++ratio) {
            U_neighbors[u][0].merge(U_neighbors[u][ratio]);
        }
    }
    return was_designated;
}

NetworKit::node KRTEdgeDesignator::designate_edge(NetworKit::node u) {
    NetworKit::node v = NetworKit::none;
    if (!U_prim.count(u)) {
        if (!U_neighbors[u].front().empty()) {
            v = *(U_neighbors[u].front().begin());
        }
        designated[u] = v;
    } else {
        for (int ratio = 0; ratio <= t; ++ratio) {
            if (!U_neighbors[u][ratio].empty()) {
                v = *(U_neighbors[u][ratio].begin());
                break;
            }
        }
        designated[u] = v;
        if (v != NetworKit::none) {
            update_rl(v);
            if (rl[v] > erl[v]) {
                update_erl(v);
            }
        }
    }
    return v;
}

void KRTEdgeDesignator::reset_if_needed(NetworKit::node v) {
    if (v == NetworKit::none || rl[v] < t) {
        return;
    }
    do {
        reset();
    } while (std::any_of(rl.begin(), rl.end(), [this](int ratio_level) {
        return ratio_level >= t;
    }));
}

void KRTEdgeDesignator::remove_edge_and_redesignate(
        NetworKit::node u, NetworKit::node v) {
    if (remove_edge(u, v)) {
        reset_if_needed(designate_edge(u));
    }
}

long double KRTEdgeDesignator::reset() {
    auto k = t;
    const long double threshold = l / (88.0 * x);
    while (k >= 4
            && get_indexed_U(k - 3).size()
                >= ratios[k - 3] * threshold * get_indexed_U(k).size()) {
        k -= 3;
    }
    auto W = get_indexed_U(k - 1);
    for (const auto &v : get_indexed_V(k - 1)) {
        while (rl[v] >= k - 1) {
            for (const auto &u : U_prim) {
                if (designated[u] == v) {
                    designated[u] = NetworKit::none;
                    break;
                }
            }
            update_rl(v);
        }
        if (rl[v] > erl[v] || rl[v] < erl[v] - 1) {
            update_erl(v);
        }
    }
    for (const auto &u : W) {
        if (designated[u] == NetworKit::none) {
            designate_edge(u);
        }
    }
    return k;
}

void KRTEdgeDesignator::initialize(const std::optional<NetworKit::Graph> &graph) {
    initialize(graph, Parameters{});
}

void KRTEdgeDesignator::initialize(
        const std::optional<NetworKit::Graph> &graph, const Parameters &parameters) {
    NetworKit::count n = graph->numberOfNodes();
    MAX_K = 2 * n, N = (n + 1) * MAX_K, M = 0;
    U = std::vector<std::vector<NetworKit::node>>(N);
    V = std::vector<std::vector<NetworKit::node>>(N);
    deg_U = std::vector<NetworKit::count>(N, 0);
    designated = std::vector<NetworKit::node>(N, NetworKit::none);
    rl = std::vector<int>(N, 0), erl = std::vector<int>(N, 0);
    initialize_parameters(parameters);

    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        for (int k = 1; k < static_cast<int>(MAX_K); k++, M++) {
            NetworKit::node left = encode_id(u, k), right = encode_id(v, k - 1);
            U[left].push_back(right);
            V[right].push_back(left);
        }
    });
    initialize_ratios();
    initialize_neighbors();
    initialize_prim();
    for (NetworKit::node i = 0; i < N; ++i) {
        designate_edge(i);
    }
}

NetworKit::node KRTEdgeDesignator::current_edge(NetworKit::node i, int k) {
    return decode_id(designated[encode_id(i, k)]);
}

void KRTEdgeDesignator::response_adversary(NetworKit::node a, int da) {
    NetworKit::node v = encode_id(a, da);
    for (const auto &u : V[v]) {
        remove_edge_and_redesignate(u, v);
    }
}

void KRTEdgeDesignator::response_adversary(NetworKit::node a, int da, NetworKit::node b, int db) {
    NetworKit::node u = encode_id(a, da), v = encode_id(b, db);
    remove_edge_and_redesignate(u, v);
}

}  // namespace Koala
