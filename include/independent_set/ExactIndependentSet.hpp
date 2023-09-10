/*
 * ExactIndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <map>

#include <independent_set/IndependentSet.hpp>

namespace Koala {

/**
 * @ingroup independent_set
 * The base class for the exact independent set algorithms.
 *
 */
class ExactIndependentSet : public IndependentSet {
 public:
    using IndependentSet::IndependentSet;

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
};

/**
 * @ingroup independentSet
 * An abstract class for recursive algorithms
 */
class RecursiveIndependentSet : public ExactIndependentSet {
 public:
    using ExactIndependentSet::ExactIndependentSet;

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
void ExactIndependentSet::restoreElements(
    std::vector<NetworKit::node>& nodes, T& edges) {
    for (auto v : nodes) {
        graph->restoreNode(v);
    }
    for (auto e : edges) {
        graph->addEdge(e.u, e.v);
    }
}

}  // namespace Koala
