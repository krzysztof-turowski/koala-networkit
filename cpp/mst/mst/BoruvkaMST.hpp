#pragma once

#include "MinimalSpanningTree.hpp"
#include "utils.hpp"
#include <optional>

namespace MST {

    class BoruvkaMST;

    class Bookkeeping;

    class FullBranchingTree {
        /**
         * Definition:
         *
         * Implementation:
         * Directed tree, where edges are directed top->down.
         * (Root has no incoming edges, leaves have no outgoing ones).
         * Created only by BoruvkaMST.
         * If Interested only in FullBranchingTree, then call `auto fbt = BoruvkaMST(G, ,true).fullBranchingTree->get()`
         **/


        friend class BoruvkaMST;

        friend class Bookkeeping;

        Graph tree;

        explicit FullBranchingTree(const Graph& G);

    public:
        [[nodiscard]]
        const Graph& getTree() const;

        [[nodiscard]]
        node getRoot() const;

        // [0, leaves) nodes are leaves.
        const count leaves;

        [[nodiscard]]
        std::optional<node> getParent(node u) const;
    };


    class BoruvkaMST final : public MinimalSpanningTree {
    public:
        std::optional<FullBranchingTree> fullBranchingTree;
        std::optional<Graph> contractedG_iAfterLimitedSteps;
        std::optional<std::vector<NodePair>> contractedG_iGetOriginalNodesFromEdgeId;

        // storeG_iContracted: if `true`, then only `2` iterations of Boruvka algorithm will be conducted
        // and the result will be stored in `G_i`. Necessary for KktMST.
        explicit BoruvkaMST(const Graph& G, bool storeOriginalGraph = false, bool storeFullBranchingTree = false,
                            bool storeG_iContracted = false);

        ~BoruvkaMST() override = default;
    };
}
