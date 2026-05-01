/*
 * ChordalGraphRecognition.hpp
 *
 *  Created on:
 *      Author:
 */

#pragma once

#include <limits>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {


/**
 * @ingroup recognition
 * Represents a sequence of vertices. If the graph is chordal,
 * this will be a Perfect Elimination Ordering.
 */
struct PerfectEliminationOrdering {
   std::vector<NetworKit::count> alpha;     // node -> position in ordering
   std::vector<NetworKit::node> alpha_inv;  // position in ordering -> node

   PerfectEliminationOrdering(NetworKit::count n = 0, NetworKit::count bound = 0);

   void set(NetworKit::node v, NetworKit::count pos);
};

/**
 * @ingroup recognition
 * Base class for chordal graph recognition algorithms.
 */
class ChordalGraphRecognition : public NetworKit::Algorithm {
public:
   enum class State {
      UNKNOWN,
      CHORDAL,
      NOT_CHORDAL
   };

   /**
    * Given an input graph, set up the chordal graph recognition.
    *
    * @param graph The input graph.
    */
   explicit ChordalGraphRecognition(const NetworKit::Graph& graph);

   /**
    * Execute the chordal graph recognition procedure.
    */
   virtual void run() = 0;

   /**
    * Return the result found by the algorithm.
    *
    * @return true if the graph is chordal, false otherwise.
    */
   virtual bool isChordal() const;

   /**
    * Return the graph type found by the algorithm.
    *
    * @return State of the graph.
    */
   virtual State getState() const;

   /**
    * Verify the result found by the algorithm.
    */
   virtual void check() const;

protected:
   NetworKit::Graph graph;
   State is_chordal;
};

/**
 * @ingroup recognition
 * The class for recognition of chordal graphs procedure from
 * Tarjan, Yannakakis, "Simple Linear-Time Algorithms to Test Chordality of Graphs,
 * Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs".
 */
class MaximumCardinalitySearchChordalGraphRecognition : public ChordalGraphRecognition {
public:
   using ChordalGraphRecognition::ChordalGraphRecognition;

   /**
    * Execute the MCS chordal graph recognition algorithm.
    */
   void run() override;

private:
   /**
    * Phase 1: Construct an ordering using MCS.
    */
   PerfectEliminationOrdering getMCSOrdering() const;

   /**
    * Phase 2: Verify if the generated ordering has zero fill-in.
    */
   bool isZeroFillIn(const PerfectEliminationOrdering& ordering) const;

   /**
    * Helper structure to maintain the array of sets in O(1) time
    * using doubly linked lists.
    */
   class SetArray {
   public:
      SetArray(const NetworKit::Graph& graph);
      void add(NetworKit::node v, NetworKit::count set);
      void remove(NetworKit::node v, NetworKit::count set);
      NetworKit::node getHead(NetworKit::count set) const;
      bool isEmpty(NetworKit::count set) const;

      static constexpr NetworKit::count INITIAL_SET = 0;

   private:
      struct VertexNode {
         NetworKit::node prev{ NetworKit::none };
         NetworKit::node next{ NetworKit::none };
      };

      std::vector<NetworKit::node> head;
      std::vector<VertexNode> vertices;
   };
};


/**
 * @ingroup recognition
 * The class for recognition of chordal graphs procedure from
 * Rose, Tarjan, Lueker, "Algorithmic Aspects of Vertex Elimination on Graphs"
 * using Lexicographic Breadth-First Search.
 */
class LexBFSChordalGraphRecognition : public ChordalGraphRecognition {
public:
   using ChordalGraphRecognition::ChordalGraphRecognition;

   /**
    * Execute the LexBFS chordal graph recognition algorithm.
    */
   void run() override;

private:
   /**
    * Phase 1: Construct an ordering using LexBFS.
    */
   PerfectEliminationOrdering getLexBFSOrdering() const;

   /**
    * Phase 2: Verify if the generated ordering has zero fill-in
    */
   bool isZeroFillIn(const PerfectEliminationOrdering& ordering) const;

   /**
    * Helper structure to maintain the Queue of Sets in O(1) time
    * using doubly linked lists of sets and doubly linked lists of vertices.
    */
   class SetQueue {
   public:
      SetQueue(const NetworKit::Graph& graph);
      NetworKit::node popVertex();
      void moveToNewSet(NetworKit::node w);
      void removeEmptySets();
      bool isProcessed(NetworKit::node w) const;

   private:
      static constexpr NetworKit::count NO_SET = 0;
      static constexpr NetworKit::count INITIAL_SET = 1;
      static constexpr NetworKit::count PROCESSED_SET =
         std::numeric_limits<NetworKit::count>::max();

      void removeVertex(NetworKit::node v, NetworKit::count set);
      void addVertex(NetworKit::node v, NetworKit::count set);
      void removeSet(NetworKit::count set);

      struct SetNode {
         NetworKit::node head{ NetworKit::none };
         NetworKit::count prev{ NO_SET };
         NetworKit::count next{ NO_SET };
         NetworKit::count split_into{ NO_SET };
      };

      struct VertexNode {
         NetworKit::count set_id{ NO_SET };
         NetworKit::node prev{ NetworKit::none };
         NetworKit::node next{ NetworKit::none };
      };

      std::vector<SetNode> sets;
      std::vector<VertexNode> vertices;

      NetworKit::count first_set;
      std::vector<NetworKit::count> free_sets;
      std::vector<NetworKit::count> active_splits;
   };
};


} /* namespace Koala */
