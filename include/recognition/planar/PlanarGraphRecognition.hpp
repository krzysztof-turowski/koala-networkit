/*
 * PlanarGraphRecognition.hpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

 #pragma once

 #include <optional>
 #include <stack>
 #include <vector>
 #include <list>
 #include <networkit/base/Algorithm.hpp>
 #include <networkit/components/ConnectedComponents.hpp>
 #include <networkit/graph/Graph.hpp>
 #include <networkit/graph/GraphTools.hpp>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <stack>
#include <vector>
#include <deque>
#include <stdexcept>
#include <numeric>
#include <iterator>
#include <utility>
#include <limits>
#include <array>

 
 namespace Koala {
 
     const int LEFT = 1;
     const int RIGHT = 2;
 
 /**
  * @ingroup recognition
  * The abscract class for recognition of planar graphs p.
  *
  */
     class PlanarGraphRecognition : public NetworKit::Algorithm {
     public:
         enum class State { NOT_PLANAR, PLANAR };
 
         /**
          * Given an input graph, set up the planar graph recognition.
          *
          * @param graph The input graph.
          */
         PlanarGraphRecognition(NetworKit::Graph &graph, bool embedding);
 
         /**
          * Execute the planar graph recognition procedure.
          */
         virtual void run() = 0;
 
         /**
          * Return the result found by the algorithm.
          *
          * @return true if the graph is planar, false otherwise.
          */
         PlanarGraphRecognition::State isPlanar() const;
 
         /**
          * Return the embbending found by the algorithm.
          *
          * @return pointer to embedding if the graph is planar, nullptr otherwise.
          */
         void getEmbedding() const;
 
     protected:
         NetworKit::Graph graph;
         State is_planar;
         bool embedding;
     };
 
 /**
  * @ingroup recognition
  * The class for recognition of perfect graphs procedure from
  * Chudnovsky, Cornuejols, Liu, Seymour, Vuskovic, "Recognizing Berge graphs".
  *
  */
     class HopcroftTarjan : public PlanarGraphRecognition {
     public:
         HopcroftTarjan(NetworKit::Graph &graph, bool embedding);
         void run() override;
 
     private:
         NetworKit::count N;
 
         void dfs(NetworKit::Graph &g, NetworKit::node v,
                  std::vector<bool> &reached);
         void dfs_in_make_biconnected_graph(NetworKit::Graph &g, NetworKit::node v,
                                            int &dfs_count,
                                            std::vector<bool> &reached,
                                            std::vector<NetworKit::node> &dfsnum,
                                            std::vector<NetworKit::node> &lowpt,
                                            std::vector<NetworKit::node> &parent);
 
         void make_biconnected_graph(NetworKit::Graph &g);
 
         void dfs_in_reorder(NetworKit::Graph &g, NetworKit::node v,
                             std::vector<NetworKit::Edge> &Add, int &dfs_count,
                             std::vector<bool> &reached,
                             std::vector<NetworKit::node> &dfsnum,
                             std::vector<NetworKit::node> &lowpt,
                             std::vector<NetworKit::node> &lowpt2,
                             std::vector<NetworKit::node> &parent);
 
         void reorder(NetworKit::Graph &g, std::vector<NetworKit::node> &parent,
                      std::vector<NetworKit::node> &dfsnum,
                      std::vector<NetworKit::Edge> &order);
 
         bool strong_planar(NetworKit::Graph &g, NetworKit::Edge e,
                            std::list<int> &Att, std::vector<int> &alpha,
                            std::vector<NetworKit::node> &dfsnum,
                            std::vector<NetworKit::node> &parent);
 
         PlanarGraphRecognition::State planar();
     };
     /*
     class BoyerMyrvold : public PlanarGraphRecognition {
     public:
       //  BoyerMyrvold(NetworKit::Graph &graph, bool embedding);
       //  void run() override;

     };
     */
 
 }

 
 namespace Koala {
    using namespace NetworKit;

    class LeftRightPlanarity {
       public:
        struct FringeOpposedSubset {
            std::deque<int> left, right;
    
            FringeOpposedSubset(int h) {
                left.push_back(h);
            }
    
            int l_lo() const { return left.back(); }
            int l_hi() const { return left.front(); }
            int r_lo() const { return right.back(); }
            int r_hi() const { return right.front(); }
    
            bool hasRight() const { return !right.empty(); }
            bool hasLeft() const { return !left.empty(); }
        };
    
        struct Fringe {
            std::deque<FringeOpposedSubset> fops;
    
            Fringe() {}
            Fringe(int dfs_h) { fops.emplace_back(dfs_h); }
    
            FringeOpposedSubset& H();
            FringeOpposedSubset& L();
    
            bool operator<(const Fringe& other) const;
            void merge(Fringe& other);
            void prune(int dfs_height);
            std::pair<bool, bool> lrCondition(int dfs_height);
            void swapSide();
            void mergeTAlikeEdges();
            void mergeTOppositeEdgesInto(Fringe& other);
            void alignDuplicates(int dfs_h);
            void makeOnionStructure(Fringe& other);
        };
    

        bool isPlanar(const Graph& G);
    private:
        bool lrAlgorithm(const Graph& G, node root, std::vector<int>& dfsHeight);

        void mergeFringes(std::vector<std::vector<Fringe>>& fringes, int dfs_height);
    
        Fringe getMergedFringe(std::vector<Fringe>& upperFringes);
    };
    
} // namespace PlanarityTester
