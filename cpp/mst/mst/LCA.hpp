#pragma once

#include "BoruvkaMST.hpp"
#include "utils.hpp"


namespace MST {
    class LCA {
        /*
         * This class implements <O(n), O(1)> algorithm from Bender & Farach-Colton "The LCA Problem Revisited" (2000).
         * Since NetworKit doesn't provide a separate class for a tree, we use FullBranchingTree. Note that the
         * algorithm in general is suitable for any tree.
         */

    public:
        using depth_t = count;
        using mask_t = index;

        // preprocessing, O(|V|)
        explicit LCA(const FullBranchingTree& fbt);

        // query, O(1)
        [[nodiscard]]
        node lca(node u, node v) const;

    protected:
        struct QueryToBlockInfos;
        static const count MIN_ALLOWED_BLOCK_SIZE = 2;
        const FullBranchingTree& fbt;

        // E[] in the paper. E[i] is the label of the ith node visited in the Euler tour of T.
        std::vector<node> eulerTour;

        // L[] in the paper. Distance from root. L[i] is the depth of node E[i].
        std::vector<depth_t> depth;

        // R in the paper; "time-in" in DFS nomenclature. E[R[i]] = i;
        // note that E[] could have many entries for `i` - R[i] points to the first one.
        std::vector<index> representative;


        count blockSize; // log(eulerTour.size())/2
        /*
        * precomputedRMQSmallBlocks[mask][l][c]: for a given `mask` (defined by +-1 sequence depicting depths using
         * mapping {-1 -> 0; 1 -> 1}) stores argMin of <sequence that generated given mask> in interval [l, l+c] (both inclusive).
        * i.e.: precomputedRMQSmallBlocks[mask][0][blockSize-1] will give an answer for the whole block.
        */
        std::vector<std::vector<std::vector<index>>> precomputedRMQSmallBlocks;

        // --- BLOCK ARRAYS --- //
        // A' in the paper. - depthMinOfBlock[i] holds minimum element of depth[i*blockSize...i*(blockSize+1)-1]
        std::vector<depth_t> depthMinOfBlock;

        // similar to B in the paper. Difference: in our implementation, eulerTourIndexRealizingDepthMinOfBlock[]
        // holds values incremented by `blockOffset`
        std::vector<index> eulerTourIndexRealizingDepthMinOfBlock;

        // depthBlock may be described as (startingValue, mask), where mask is a +-1s sequence.
        // The vector below provides mapping: depth[i*blockSize...i*(blockSize+1)-1] -> mask
        std::vector<mask_t> depthBlockToMask;

        // --- SPARSE TABLE --- //
        // M[,] in the paper. sparseTableDepthMinOfBlockRMQ[l][j] is the index of minimum element in depthMinOfBlock  i [l...l+2^j)
        std::vector<std::vector<index>> sparseTableDepthMinOfBlockRMQ;

        void eulerTourDFS(node u, count level);

        void populatePrecomputedRMQSmallBlocks();

        void populateBlockArrays();

        void populateSparseTable();

        [[nodiscard]]
        static mask_t getMask(std::vector<depth_t>::iterator begin, std::vector<depth_t>::iterator end);

        [[nodiscard]]
        index querySmallBlock(mask_t mask, index l, index r) const;

        [[nodiscard]]
        node lcaNaive(node u, node v) const;
    };
}
