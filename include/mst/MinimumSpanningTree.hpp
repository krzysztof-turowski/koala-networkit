/*
 * MinimumSpanningTree.hpp
 *
 *  Created on: 08.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace Koala {

/**
 * @ingroup mst
 * The base class for the minimum spanning tree algorithms.
 *
 */
class MinimumSpanningTree : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the minimum spanning tree procedure.
     *
     * @param graph The input graph.
     */
    explicit MinimumSpanningTree(NetworKit::Graph &graph);

    /**
     * Return the spanning tree found by the algorithm.
     *
     * @return a spanning tree.
     */
    const NetworKit::Graph& getForest() const;

    /**
     * Verify the result found by the algorithm using O(n + m) MST verification algorithm
     * from Hagerup, An Even Simpler Linear-Time Algorithm for Verifying Minimum Spanning Trees.
     */
    void check() const;

 protected:
    using NodePair = std::pair<NetworKit::node, NetworKit::node>;

    std::optional<NetworKit::Graph> graph, tree;
};

/**
 * @ingroup mst
 * The class for the Kruskal minimum spanning tree algorithm
 */
class KruskalMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Kruskal minimum spanning tree algorithm.
     */
    void run();
};

/**
 * @ingroup mst
 * The class for the Prim minimum spanning tree algorithm
 */
class PrimMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Prim minimum spanning tree algorithm.
     */
    void run();
};

/**
 * @ingroup mst
 * The class for the Boruvka minimum spanning tree algorithm
 */
class BoruvkaMinimumSpanningTree : public MinimumSpanningTree {
    friend class MinimumSpanningTree;
    friend class Chazelle2000MinimumSpanningTree;
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Boruvka minimum spanning tree algorithm.
     */
    void run();

 protected:
    static std::optional<NetworKit::Graph> iterate(
        NetworKit::Graph &G, NetworKit::Graph &F,
        NetworKit::UnionFind &union_find, std::map<NodePair, NodePair> &E,
        NetworKit::count steps, bool get_branching_tree);
};

/**
 * @ingroup mst
 * The class for the Karger-Klein-Tarjan randomized minimum spanning tree algorithm
 */
class KargerKleinTarjanMinimumSpanningTree final : public BoruvkaMinimumSpanningTree {
 public:
    using BoruvkaMinimumSpanningTree::BoruvkaMinimumSpanningTree;

    /**
     * Execute the Karger-Klein-Tarjan randomized minimum spanning tree algorithm.
     */
    void run();

 protected:
    static void recurse(NetworKit::Graph &G, NetworKit::Graph &F);
    static void discard_random_edges(NetworKit::Graph &G, NetworKit::Graph &subgraph);
    static void remove_heavy_edges(NetworKit::Graph &G, NetworKit::Graph &subgraph);
};

class ChazelleRubinfeldTrevisanMinimumSpanningTree final : public MinimumSpanningTree{
public:
   using MinimumSpanningTree::MinimumSpanningTree;
   /**
    * Execute Chazelle-Rubinfeld-Trevisan randomized minimum spanning tree weight algorithm.
    * 
    * eps - a constant in (0, 0.5), used to bound running time and result accuracy.
    * The smaller the eps, the more accurate the output.
    * 
    * w - maximum edge weight - the algorithm assumes all edges have weights from {1,...,w}
    * as the algorithm is sublinear, w cannot be determined at runtime.
    */
   void run(unsigned int w, float eps = 0.1);

   /**
    * Do not use this function, instead use run with parameters
    * This function is neccessary for the class to compile, but is not implemented.
    */
   void run();

   /**
    * The algorithm does not calculate minimum spanning tree, only it's approximate weight.
    * Thus getForest method throws an exception.
    */
   const NetworKit::Graph& getForest() const;
   
   /** 
    * Get the approximate weight of minimum spanning tree, calculated in run() method.
   */
   float getTreeWeight() const; 
private:
   float calculateApproximateDegree(float eps) const;
   float calculateApproximateCCsCount(float eps, int bfs_bound, unsigned int w, unsigned int w_bound) const; 
   float calculateApproximateTreeWeight(float eps, unsigned int w) const;
   float treeWeight = -1;
};

class Chazelle2000MinimumSpanningTree final : public MinimumSpanningTree {
public:
   using MinimumSpanningTree::MinimumSpanningTree;
   void run();

private:
   NetworKit::Graph msf(NetworKit::Graph G, std::map<NodePair, int>& edgeId, int t);
   std::tuple<NetworKit::Graph, const std::map<NodePair, int>&, NetworKit::Graph> boruvkaSteps(NetworKit::Graph G, std::map<NodePair, int>& edgeId, int c);

   struct EdgeInfo {
      bool bad = false;
      bool in_msf = false;
      int original_cost = 0;
   };
   std::vector<EdgeInfo> edges; 
   static constexpr int MIN_NUMBER_NODES = 10;
   static constexpr int C = 10;
};

}  /* namespace Koala */
