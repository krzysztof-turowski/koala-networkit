/*
 * KrtEdgeDesignator.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <map>
#include <optional>
#include <unordered_set>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

class KRTEdgeDesignator {
    // Constraint: R0 * L/X >= 176, then T = ceil(log(N) / log((R0 * L)/(88 * X))) + 4;
    // Parameters selected empirically to (1) match constraints, (2) get best results on datasets
    static constexpr long double R0 = 0.7, X = 2;
    static constexpr int L = 512;
    static constexpr int T = 7;

    int N, M, MAX_K;
    std::vector<int> degU, designated, rl, erl;
    std::vector<long double> ratios;
    std::unordered_set<int> U_prim, V_prim;
    std::vector<std::vector<int>> U, V;
    std::vector<std::vector<std::unordered_set<int>>> U_neighbors;

    void initialize_prim();
    void initialize_ratios();
    void initialize_neighbors();

    std::unordered_set<int> get_indexed_U(int);
    std::unordered_set<int> get_indexed_V(int);

    int encodeId(NetworKit::node, int) const;
    NetworKit::node decodeId(int) const;

    void update_rl(int);
    void update_erl(int);

    void remove_edge(int, int);
    int designate_edge(int);

    long double reset();

 public:
    void initialize(const std::optional<NetworKit::Graph>&);
    NetworKit::node current_edge(NetworKit::node, int);
    void response_adversary(NetworKit::node, int);
    void response_adversary(NetworKit::node, int, NetworKit::node, int);
};

} /* namespace Koala */
