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

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * Residual-edge designator used by KingRaoTarjanMaximumFlow.
 *
 * The data structure tracks one candidate residual edge per encoded vertex-level state and
 * updates designations in response to relabels and rejected edges.
 */
class KRTEdgeDesignator {
 public:
    /**
     * Parameters controlling ratio levels and reset thresholds.
     */
    struct Parameters {
        // The general strategy does not prescribe R0 or X. Defaults are 0.7 and 2.
        // If omitted, L is the smallest value satisfying R0 * L / X > 176 and T is
        // computed from the paper's formula.
        std::optional<long double> R0 = std::nullopt;
        std::optional<NetworKit::count> L = std::nullopt;
        std::optional<long double> X = std::nullopt;
        std::optional<int> T = std::nullopt;
    };

 private:
    long double r0, x;
    NetworKit::count l;
    int t;

    NetworKit::count N, M, MAX_K;
    std::vector<NetworKit::count> deg_U;
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

    NetworKit::node encode_id(NetworKit::node, int) const;
    NetworKit::node decode_id(NetworKit::node) const;

    void update_rl(NetworKit::node);
    void update_erl(NetworKit::node);

    bool remove_edge(NetworKit::node, NetworKit::node);
    NetworKit::node designate_edge(NetworKit::node);
    void reset_if_needed(NetworKit::node);
    void remove_edge_and_redesignate(NetworKit::node, NetworKit::node);

    long double reset();

 public:
    /**
     * Initialize the designator with default parameters.
     *
     * @param graph Residual graph represented as an optional graph.
     */
    void initialize(const std::optional<NetworKit::Graph>&);

    /**
     * Initialize the designator with explicit parameters.
     *
     * @param graph Residual graph represented as an optional graph.
     * @param parameters Designator tuning parameters.
     */
    void initialize(const std::optional<NetworKit::Graph>&, const Parameters&);

    /**
     * Return the currently designated edge endpoint for a vertex-level state.
     *
     * @return The designated endpoint, or NetworKit::none if no edge is designated.
     */
    NetworKit::node current_edge(NetworKit::node, int);

    /**
     * Notify the designator that a vertex has been relabeled.
     */
    void response_adversary(NetworKit::node, int);

    /**
     * Notify the designator that an edge candidate was rejected.
     */
    void response_adversary(NetworKit::node, int, NetworKit::node, int);
};

}  // namespace Koala
