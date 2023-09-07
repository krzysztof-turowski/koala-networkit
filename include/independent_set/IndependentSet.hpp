/*
 * IndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <map>
#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup independent_set
 * The base class for the independent set problem algorithms.
 *
 */
class IndependentSet : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the independent set problem procedure.
     *
     * @param graph The input graph.
     */
    explicit IndependentSet(const NetworKit::Graph &graph);

    /**
     * Return the independent set found by the algorithm.
     *
     * @return a set of nodes that form an independent set.
     */
    const std::set<NetworKit::node>& getIndependentSet() const;

    /**
     * Execute the maximum independent set finding procedure.
     */
    virtual void run() = 0;

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

 protected:
    /**
     * Comparator for EdgeSet type
     */
    static bool edgeComparator(const NetworKit::Edge& a, const NetworKit::Edge& b);
    using EdgeSet = std::set<NetworKit::Edge, decltype(&edgeComparator)>;

    /**
     * @return N(v) - neigbors of vertex v
    */
    std::vector<NetworKit::node> getNeighbors(NetworKit::node v) const;

    /**
     * @return N[v] - neigbors of vertex v and v
    */
    std::vector<NetworKit::node> getNeighborsPlus(NetworKit::node v) const;

    /**
     * @return N^2(v) - neigbors in a distance 2 of v
    */
    std::vector<NetworKit::node> getNeighbors2(NetworKit::node v) const;

    EdgeSet getConnectedEdges(std::vector<NetworKit::node>& nodes) const;
    EdgeSet getInducedEdges(std::vector<NetworKit::node>& nodes) const;
    std::vector<NetworKit::node> getMirrors(NetworKit::node v) const;

    NetworKit::node getMinimumDegreeNode() const;
    NetworKit::node getMaximumDegreeNode() const;
    void removeElements(std::vector<NetworKit::node> nodes);
    template <typename T>
    void restoreElements(std::vector<NetworKit::node>& nodes, T& edges);
    std::vector<NetworKit::node> runIndependentSetDegree2() const;

    std::optional<NetworKit::Graph> graph;
    std::set<NetworKit::node> independentSet;
};

/**
 * @ingroup independentSet
 * The class for the brute force algorithm.
 */
class BruteForceIndependentSet final : public IndependentSet {
 public:
    using IndependentSet::IndependentSet;

    /**
     * Execute the brute force algorithm.
     */
    void run();
};

/**
 * @ingroup independentSet
 * An abstract class for recursive algorithms
 */
class RecursiveIndependentSet : public IndependentSet {
 public:
    using IndependentSet::IndependentSet;

    /**
     * Start recursive calls and extract result
     */
    virtual void run();

    /**
     * Actual recursive function
     */
    virtual std::vector<NetworKit::node> recursive() = 0;
};

/**
 * @ingroup independentSet
 * The class for the mis1 algorithm.
 */
class Mis1IndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the mis1 algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis2 algorithm.
 */
class Mis2IndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the mis2 algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis3 algorithm.
 */
class Mis3IndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the mis3 algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis4 algorithm.
 */
class Mis4IndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the mis4 algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis5 algorithm.
 */
class Mis5IndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the mis5 algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the Measue and Conquer Simple O(2^0.288n) algorithm.
 */
class MeasureAndConquerIndependentSet final : public RecursiveIndependentSet {
 public:
    using RecursiveIndependentSet::RecursiveIndependentSet;

 private:
    /**
     * Execute the Measue and Conquer Simple algorithm.
     */
    std::vector<NetworKit::node> recursive();
};

template <typename T>
void IndependentSet::restoreElements(
    std::vector<NetworKit::node>& nodes, T& edges) {
    for (auto v : nodes) {
        graph->restoreNode(v);
    }
    for (auto e : edges) {
        graph->addEdge(e.u, e.v);
    }
}

} /* namespace Koala */
