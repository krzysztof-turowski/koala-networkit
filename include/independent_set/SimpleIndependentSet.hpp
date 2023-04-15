/*
 * SimpleIndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <functional>
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
class SimpleIndependentSet : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the independent set problem procedure.
     *
     * @param graph The input graph.
     */
    SimpleIndependentSet(const NetworKit::Graph &graph);

    /**
     * Return the independent set found by the algorithm.
     *
     * @return a map from nodes to (true <=> belongs to the independent set).
     */
    const std::map<NetworKit::node, bool>& getIndependentSet() const;

protected:
    std::function<bool(NetworKit::Edge a, NetworKit::Edge b)> edgeComparator;
    using EdgeSet = std::set<NetworKit::Edge, decltype(edgeComparator)>;

    std::vector<NetworKit::node> getNeighbors(NetworKit::node v) const;
    std::vector<NetworKit::node> getNeighborsPlus(NetworKit::node v) const;
    std::set<NetworKit::node> getNeighbors2(NetworKit::node v) const;
    std::set<NetworKit::node> getNeighbors2Plus(NetworKit::node v) const;
    EdgeSet getConnectedEdges(std::vector<NetworKit::node>& nodes) const;
    EdgeSet getInducedEdges(std::vector<NetworKit::node>& nodes) const;
    EdgeSet getAllEdges() const;
    std::vector<NetworKit::node> getAllNodes() const;
    std::vector<NetworKit::node> getMirrors(NetworKit::node v) const;
    void dfs(NetworKit::node v, std::vector<bool>& visited);

    NetworKit::count getGraphsMaximumDegree() const;
    NetworKit::node getMinimumDegreeNode() const;
    NetworKit::node getMaximumDegreeNode() const;
    void removeElements(std::vector<NetworKit::node> nodes);
    void restoreElements(
        std::vector<NetworKit::node>& nodes, 
        EdgeSet& edges);
    std::vector<NetworKit::node> runIndependentSetDegree2() const;

    const std::optional<NetworKit::Graph> graph;
    std::map<NetworKit::node, bool> independentSet;
};

/**
 * @ingroup independentSet
 * The class for the brute force algorithm.
 */
class BruteForceIndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the brute force algorithm.
     */
    void run();
};

/**
 * @ingroup independentSet
 * The class for the mis1 algorithm.
 */
class Mis1IndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;
    
    /**
     * Execute the mis1 algorithm.
     */
    void run();

private:
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis2 algorithm.
 */
class Mis2IndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the mis2 algorithm.
     */
    void run();
private:
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis3 algorithm.
 */
class Mis3IndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the mis3 algorithm.
     */
    void run();
private:
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis4 algorithm.
 */
class Mis4IndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the mis4 algorithm.
     */
    void run();
private:
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the mis5 algorithm.
 */
class Mis5IndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the mis5 algorithm.
     */
    void run();
private:
    std::vector<NetworKit::node> recursive();
};

/**
 * @ingroup independentSet
 * The class for the Measue and Conquer Simple O(2^0.288n) algorithm.
 */
class MeasureAndConquerIndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the mis5 algorithm.
     */
    void run();
private:
    std::vector<NetworKit::node> recursive();
};


} /* namespace Koala */
