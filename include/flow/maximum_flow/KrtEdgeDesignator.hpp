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
 public:
    struct Parameters {
        // The general strategy does not prescribe r0 or x. Defaults are 0.7 and 2.
        // If omitted, l is the smallest value satisfying r0 * l / x > 176 and t is
        // computed from the paper's formula.
        std::optional<long double> r0 = std::nullopt;
        std::optional<NetworKit::count> l = std::nullopt;
        std::optional<long double> x = std::nullopt;
        std::optional<int> t = std::nullopt;
    };

 private:
    long double r0, x;
    NetworKit::count l;
    int t;

    NetworKit::count N, M, MAX_K;
    std::vector<NetworKit::count> degU;
    std::vector<NetworKit::node> designated;
    std::vector<int> rl, erl;
    std::vector<long double> ratios;
    std::unordered_set<NetworKit::node> U_prim, V_prim;
    std::vector<std::vector<NetworKit::node>> U, V;
    std::vector<std::vector<std::unordered_set<NetworKit::node>>> U_neighbors;

    void initialize_prim();
    void initialize_parameters(const Parameters&);
    void initialize_ratios();
    void initialize_neighbors();

    std::unordered_set<NetworKit::node> get_indexed_U(int);
    std::unordered_set<NetworKit::node> get_indexed_V(int);

    NetworKit::node encodeId(NetworKit::node, int) const;
    NetworKit::node decodeId(NetworKit::node) const;

    void update_rl(NetworKit::node);
    void update_erl(NetworKit::node);

    bool remove_edge(NetworKit::node, NetworKit::node);
    NetworKit::node designate_edge(NetworKit::node);
    void reset_if_needed(NetworKit::node);
    void remove_edge_and_redesignate(NetworKit::node, NetworKit::node);

    long double reset();

 public:
    void initialize(const std::optional<NetworKit::Graph>&);
    void initialize(const std::optional<NetworKit::Graph>&, const Parameters&);
    NetworKit::node current_edge(NetworKit::node, int);
    void response_adversary(NetworKit::node, int);
    void response_adversary(NetworKit::node, int, NetworKit::node, int);
};

} /* namespace Koala */
