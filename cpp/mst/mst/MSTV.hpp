#pragma once

#include "MinimalSpanningTree.hpp"
#include "BoruvkaMST.hpp"

namespace MST {
    /*
     * O(V+E) Minimal Spanning Tree Verification Algorithm. Returns true if a given mst is truly minimal.
     * Torben Hagerup, 2010: "An Even Simpler Linear-Time Algorithm for Verifying Minimum Spanning Trees"
     */
    bool HagerupMSTV(const MinimalSpanningTree& mst);

    /*
     * Similar to HagerupMSTV.
     * F is a spanning tree of G.
     * For each edge e=(u,v) in G and a path P=<u, ..., v> in F, this function removes e if w(e) > w(e') for e' in P.
     */
    void removeF_HeavyEdges(Graph& G, const Graph& F);
}
