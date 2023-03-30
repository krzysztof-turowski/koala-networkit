/*
 * MaximumFlow.hpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <map>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <flow/maximum_flow/DynamicTree.hpp>
#include <flow/maximum_flow/KrtEdgeDesignator.hpp>

#define DBG if(1)

class Logger {
    std::string tag = "KRT";
public:

    void log(std::string msg, std::pair<int, int> e, int val) {
        return log(msg, e.first, e.second, val);
    }

    void log(std::string msg, int u, int v, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " " + "(" << u << ", " << v << ") = " << val << std::endl;
    }

    void log(std::string msg, int u, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " " + "(" << u << ") = " << val << std::endl;
    }

    void log(std::string msg, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " = " << val << std::endl;
    }

    void log(std::string msg) {
        DBG std::cerr << "[" + tag + "]" + " " + msg << std::endl;
    }
};

namespace Koala {

/**
 * @ingroup flow
 * The base class for the max flow algorithms.
 *
 */
class MaximumFlow : public NetworKit::Algorithm {
public:
    /**
     * Given an input graph, set up the greedy vertex coloring procedure.
     *
     * @param graph The input graph.
     * @param s     The source vertex.
     * @param t     The sink vertex.
     */
    MaximumFlow(const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t);

    /**
     * Return the flow size found by the algorithm.
     *
     * @return a total flow value.
     */
    int getFlowSize() const;

protected:
    const std::optional<NetworKit::Graph> graph;
    NetworKit::node source, target;
    std::map<std::pair<NetworKit::node, NetworKit::node>, int> flow;
    int flow_size;
    
    static std::pair<NetworKit::node, NetworKit::node> reverse(
        const std::pair<NetworKit::node, NetworKit::node>&);
    
    virtual void initialize() = 0;

    Logger logger;
};

/**
 * @ingroup flow
 * The class for the King-Rao-Tarjan maximum flow algorithm
 */
class KingRaoTarjanMaximumFlow final : public MaximumFlow {

public:
    using MaximumFlow::MaximumFlow;

    /**
     * Execute the King-Rao-Tarjan maximum flow algorithm.
     */
    void run();

private:
    std::map<std::pair<NetworKit::node, NetworKit::node>, int> capacity;
    std::map<NetworKit::node, int> d, excess, hidden_excess;
    std::set<std::pair<NetworKit::node, NetworKit::node>> E_star;
    std::set<int> positive_excess;

    DynamicTree dynamic_tree;
    KRTEdgeDesignator edge_designator;

    int get_visible_excess(NetworKit::node);
    void update_positive_excess_map(NetworKit::node);
    NetworKit::node get_positive_excess_node();
    int get_flow(const std::pair<NetworKit::node, NetworKit::node>&);
    void set_flow(const std::pair<NetworKit::node, NetworKit::node>&, int);
    void saturate(const std::pair<NetworKit::node, NetworKit::node>&);
    void initialize();
    std::vector<std::pair<NetworKit::node, NetworKit::node>> get_edges_list();
    void cut(const std::pair<NetworKit::node, NetworKit::node>&);
    void tree_push(NetworKit::node, NetworKit::node);
    void relabel(NetworKit::node);
    void add_edge(const std::pair<NetworKit::node, NetworKit::node>&);
};

} /* namespace Koala */
