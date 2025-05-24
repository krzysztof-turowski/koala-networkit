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
#include <queue>
#include <unordered_map> 
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <flow/maximum_flow/DynamicTree.hpp>
#include <flow/maximum_flow/KrtEdgeDesignator.hpp>

namespace Koala {

/**
 * @ingroup flow
 * The base class for the max flow algorithms.
 *
 */

 struct pair_hash {
   template <class T1, class T2>
   std::size_t operator()(const std::pair<T1, T2>& p) const {
       std::size_t h1 = std::hash<T1>{}(p.first);
       std::size_t h2 = std::hash<T2>{}(p.second);
       return h1 ^ (h2 << 1); // or use boost::hash_combine
   }
};


class MaximumFlow : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the greedy vertex coloring procedure.
     *
     * @param graph The input graph.
     * @param s     The source vertex.
     * @param t     The sink vertex.
     */
    MaximumFlow(NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t);

    /**
     * Return the flow size found by the algorithm.
     *
     * @return a total flow value.
     */
    int getFlowSize() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    NetworKit::node source, target;
    std::unordered_map<std::pair<NetworKit::node, NetworKit::node>, int,pair_hash> flow;
    int flow_size;
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
    std::set<int> positive_excess;
    std::set<std::pair<NetworKit::node, NetworKit::node>> E_star;

    DynamicTree dynamic_tree;
    KRTEdgeDesignator edge_designator;

    int get_visible_excess(NetworKit::node);
    NetworKit::node get_positive_excess_node();
    void update_positive_excess(NetworKit::node);

    int get_flow(const std::pair<NetworKit::node, NetworKit::node>&);
    void set_flow(const std::pair<NetworKit::node, NetworKit::node>&, int);
    void saturate(const std::pair<NetworKit::node, NetworKit::node>&);
    void add_edge(const std::pair<NetworKit::node, NetworKit::node>&);
    void cut(const std::pair<NetworKit::node, NetworKit::node>&);

    void initialize();
    std::vector<std::pair<NetworKit::node, NetworKit::node>> get_edges_list();
    void tree_push(NetworKit::node, NetworKit::node);
    void relabel(NetworKit::node);
};

class PushRelabel final : public MaximumFlow {
public: 
   using MaximumFlow::MaximumFlow;
   void run();
private:
   int V;
   std::unordered_map<NetworKit::node,int> height,nextedge,excess;
   std::queue<NetworKit::node> q;
   std::unordered_map<std::pair<NetworKit::node, NetworKit::node>, int, pair_hash> capacity;

   std::pair<NetworKit::node, NetworKit::node> rev(const std::pair<NetworKit::node, NetworKit::node> &);

   void push(const NetworKit::node&, const std::pair<NetworKit::node, NetworKit::node>&);
   void relabel(const NetworKit::node&);
   void discharge(const NetworKit::node&);
   void initialize();

};

class MKMFlow final : public MaximumFlow {
public: 
   using MaximumFlow::MaximumFlow;
   void run();
private:
   int V;
   std::unordered_map<NetworKit::node,int> level;
   std::unordered_map<NetworKit::node,int> inPotential, outPotential;
   std::unordered_map<std::pair<NetworKit::node, NetworKit::node>, int, pair_hash> capacity;

   std::pair<NetworKit::node, NetworKit::node> rev(const std::pair<NetworKit::node, NetworKit::node> &);

   bool buildLevelGraph(); 
   void computePotential();
   void pushForward(NetworKit::node, int);
   void pushBackward(NetworKit::node, int);
   void deleteNode(NetworKit::node);
   void initialize();
};

class BKFlow final : public MaximumFlow {
public: 
   using MaximumFlow::MaximumFlow;
   void run();
private:
   int V;
   NetworKit::node spath,tpath;
   std::unordered_map<NetworKit::node,NetworKit::node> parent;
   std::unordered_map<NetworKit::node,int> tree;
   std::unordered_map<std::pair<NetworKit::node, NetworKit::node>, int, pair_hash> capacity;
   std::queue<NetworKit::node> active;
   std::queue<NetworKit::node> orphan;

   std::pair<NetworKit::node, NetworKit::node> rev(const std::pair<NetworKit::node, NetworKit::node> &);

   int tree_capacity(NetworKit::node, NetworKit::node);
   void initialize();
   bool grow();
   int augment();
   void adopt();
   bool origin(NetworKit::node);
};

}  /* namespace Koala */
