/*
 * ExactVertexColoring.hpp
 *
 * Created on: 06.02.2023
 *   Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
 */

#pragma once

#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

class EnumerationVertexColoring : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the enumeration vertex coloring procedure.
     *
     * @param graph The input graph.
     */
    explicit EnumerationVertexColoring(const NetworKit::Graph& graph);

    /**
     * Return the coloring found by the algorithm.
     *
     * @return a map from nodes to colors.
     */
    const std::map<NetworKit::node, int> getColoring() const;

 protected:
    const std::optional<NetworKit::Graph> graph;
    std::vector<NetworKit::node> ordering;
    std::unordered_map<NetworKit::node, int> position;
    std::vector<int> current_solution, best_solution;
    std::vector<std::set<int>> feasible_colors;
    std::set<int, std::greater<int>> current_predecessors;
    int lower_bound, upper_bound, current_bound;
    int r;

    void forwards();
    void backwards();
    void determine_feasible_colors(int i);
    virtual void determine_current_predecessors(int r) = 0;
};

class BrownEnumerationVertexColoring : public EnumerationVertexColoring {
 public:
    using EnumerationVertexColoring::EnumerationVertexColoring;

    void run();

 protected:
    std::vector<NetworKit::node> greedy_largest_first_ordering();
    void determine_current_predecessors(int r);
};

class ChristofidesEnumerationVertexColoring : public BrownEnumerationVertexColoring {
 public:
    using BrownEnumerationVertexColoring::BrownEnumerationVertexColoring;

    void run();
    const std::vector<std::vector<bool>>& getTransitiveClosure() const;

 protected:
    std::vector<std::vector<bool>> transitive_closure;

    void calculate_transitive_closure();
    void determine_current_predecessors(int r);
};

class BrelazEnumerationVertexColoring : public EnumerationVertexColoring {
 public:
    using EnumerationVertexColoring::EnumerationVertexColoring;

    void run();

 protected:
    std::vector<NetworKit::node> saturation_largest_first_with_interchange();
    bool is_interchangeable(
        std::vector<int>& K, NetworKit::node new_node, std::map<NetworKit::node, int>& solution);
    std::vector<NetworKit::node> interchange_component(
        std::vector<NetworKit::node>& subgraph, std::map<NetworKit::node, int>& solution,
        NetworKit::node new_node, int alpha);
    std::vector<int> get_representatives_of_adjacent_predecessors(int i);
    void determine_current_predecessors(int r);
    void backwards();
};

class KormanEnumerationVertexColoring : public BrownEnumerationVertexColoring {
 public:
    using BrownEnumerationVertexColoring::BrownEnumerationVertexColoring;

    void run();

 protected:
    std::vector<int> new_ordering;

    void dynamic_rearrangement(int i);
    void forwards();
    void backwards();
    void determine_feasible_colors(int i, std::unordered_set<int> blocked_colors);
};

} /* namespace Koala */
