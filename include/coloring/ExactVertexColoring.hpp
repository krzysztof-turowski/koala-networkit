/*
 * ExactVertexColoring.hpp
 *
 * Created on: 06.02.2023
 *   Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
 */

#pragma once

#include <functional>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <coloring/VertexColoring.hpp>

namespace Koala {

class EnumerationVertexColoring : public VertexColoring {
 public:
    using VertexColoring::VertexColoring;

 protected:
    std::vector<NetworKit::node> ordering;
    std::unordered_map<NetworKit::node, NetworKit::count> position;
    std::vector<NetworKit::count> current_solution, best_solution;
    std::vector<std::set<NetworKit::count>> feasible_colors;
    std::set<NetworKit::count, std::greater<NetworKit::count>> current_predecessors;
    NetworKit::count lower_bound, upper_bound, current_bound;
    NetworKit::count r;

    void forwards();
    void backwards();
    void determine_feasible_colors(NetworKit::count i);
    void store_coloring();
    virtual void determine_current_predecessors(NetworKit::count r) = 0;
};

class BrownEnumerationVertexColoring : public EnumerationVertexColoring {
 public:
    using EnumerationVertexColoring::EnumerationVertexColoring;

    void run();

 protected:
    std::vector<NetworKit::node> greedy_largest_first_ordering();
    void determine_current_predecessors(NetworKit::count r);
};

class ChristofidesEnumerationVertexColoring : public BrownEnumerationVertexColoring {
 public:
    using BrownEnumerationVertexColoring::BrownEnumerationVertexColoring;

    void run();
    const std::vector<std::vector<bool>>& getTransitiveClosure() const;

 protected:
    std::vector<std::vector<bool>> transitive_closure;

    void calculate_transitive_closure();
    void determine_current_predecessors(NetworKit::count r);
};

class BrelazEnumerationVertexColoring : public EnumerationVertexColoring {
 public:
    using EnumerationVertexColoring::EnumerationVertexColoring;

    void run();

 protected:
    std::vector<NetworKit::node> saturation_largest_first_with_interchange();
    bool is_interchangeable(
        std::vector<NetworKit::count>& K, NetworKit::node new_node,
        std::map<NetworKit::node, NetworKit::count>& solution);
    std::vector<NetworKit::node> interchange_component(
        std::vector<NetworKit::node>& subgraph,
        std::map<NetworKit::node, NetworKit::count>& solution,
        NetworKit::node new_node, NetworKit::count alpha);
    std::vector<NetworKit::count> get_representatives_of_adjacent_predecessors(
        NetworKit::count i);
    void determine_current_predecessors(NetworKit::count r);
    void backwards();
};

class KormanEnumerationVertexColoring : public BrownEnumerationVertexColoring {
 public:
    using BrownEnumerationVertexColoring::BrownEnumerationVertexColoring;

    void run();

 protected:
    std::vector<NetworKit::count> new_ordering;

    void dynamic_rearrangement(NetworKit::count i);
    void forwards();
    void backwards();
    void determine_feasible_colors(
        NetworKit::count i, std::unordered_set<NetworKit::count> blocked_colors);
};

} /* namespace Koala */
