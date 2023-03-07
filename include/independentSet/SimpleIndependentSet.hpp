/*
 * SimpleIndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

class LightWeightGraph
{
public:
    LightWeightGraph(const NetworKit::Graph &graph);
    void hide(int v);
    void unhide(int v);
    void hide (std::vector<int>& v);
    void unhide (std::vector<int>& v);
    std::vector<int> n(int v);
    int lowestDegVerticle(); // assume at least 1 size
    bool isEmpty();
    void print();
    
private:
    std::vector<std::vector<int>> adj;
    std::vector<bool> hidden;
};

/**
 * @ingroup independentSet
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
    const std::optional<NetworKit::Graph> graph;
    std::map<NetworKit::node, bool> independentSet;
    LightWeightGraph lightGraph;
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
    std::vector<int> recursive();
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
};

} /* namespace Koala */
