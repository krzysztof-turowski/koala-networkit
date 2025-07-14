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

struct EdgeHash {
    std::size_t operator()(const NetworKit::Edge& e) const {
        std::size_t h1 = std::hash<NetworKit::node>{}(e.u);
        std::size_t h2 = std::hash<NetworKit::node>{}(e.v);
        return h1 ^ (h2 << 1);  // combine hashes
    }
};

struct EdgeEqual {
    bool operator()(const NetworKit::Edge& lhs, const NetworKit::Edge& rhs) const {
        return lhs.u == rhs.u && lhs.v == rhs.v;
    }
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
      std::size_t h1 = std::hash<T1>{}(p.first);
      std::size_t h2 = std::hash<T2>{}(p.second);
      return h1 ^ (h2 << 1);
    }
};

class MaximumFlow : public NetworKit::Algorithm {
 public:
    /**
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
    int flow_size;
};

}  /* namespace Koala */
