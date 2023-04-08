#include "LCA.hpp"
#include "BoruvkaMST.hpp"
#include <bitset>
#include "utils.hpp"


namespace MST {
    struct LCA::QueryToBlockInfos {
        struct BlockInfo {
            index blockNumber;
            index indexWithinBlock;
        };
        BlockInfo left;
        BlockInfo right;

        QueryToBlockInfos(const LCA *lca, node u, node v) {
            index reprL, reprR;
            std::tie(reprL, reprR) = std::minmax(lca->representative[u], lca->representative[v]);
            left = {reprL / lca->blockSize, reprL % lca->blockSize};
            right = {reprR / lca->blockSize, reprR % lca->blockSize};
        }

        [[nodiscard]]
        bool theSameBlock() const { return left.blockNumber == right.blockNumber; }

        [[nodiscard]]
        bool adjacentBlocks() const { return right.blockNumber - left.blockNumber == 1; }
    };

    void LCA::eulerTourDFS(const node u, const depth_t level) {
        // we cannot use NetworKit's DFS because it accepts just one handle that is called during traversing adjacency
        // list before following an edge - but we need the second handle to call after coming back from the recursion.
        if (representative[u] == NetworKit::none) {
            representative[u] = eulerTour.size();
        }
        eulerTour.push_back(u);
        depth.push_back(level);
        fbt.getTree().forEdgesOf(u, [this, level](const node u, const node child) {
            eulerTourDFS(child, level + 1);
            eulerTour.push_back(u);
            depth.push_back(level);
        });
    }

//    Mask: traverse input from the left, but construct int representing mask from the right.
//    In mask, LeastSignificantBit tells us about the first difference in a sequence.
    LCA::mask_t LCA::getMask(std::vector<LCA::depth_t>::iterator begin, std::vector<LCA::depth_t>::iterator end) {
        LCA::mask_t mask = 0;
        advance(begin, 1);
        for (index i = 0; begin != end; i++, advance(begin, 1)) {
            mask |= ((*begin) > (*prev(begin))) << i;
        }
        return mask;
    }

    index LCA::querySmallBlock(index blockNumber, index l, index r) const {
        auto mask = depthBlockToMask[blockNumber];
        return precomputedRMQSmallBlocks[mask][l][r - l] + (blockNumber * blockSize);
    }

    void LCA::populatePrecomputedRMQSmallBlocks() {
        const count differentMasksCount = 1 << (blockSize - 1); // `blockSize-1` because we count +-1s between elements.
        precomputedRMQSmallBlocks.resize(differentMasksCount, std::vector<std::vector<index>>(blockSize));
        for (count mask = 0; mask < differentMasksCount; mask++) { // O(sqrt(|V|)) iterations
            for (index l = 0; l < blockSize; l++) { // O(log(|V|)) iterations
                precomputedRMQSmallBlocks[mask][l].push_back(l);
                int minVal = 0;
                index indexMinVal = l;
                int currVal = 0;
                // c - how many +-1 are considered. corresponds to mask[l...l+c]
                for (size_t c = 1; l + c < blockSize; c++) { // O(log(|V|))) iterations
                    // skip first `l` elements and include `c` +-1s looking to the right.
                    if (mask & (1 << (l + c - 1))) {
                        currVal++;
                    } else {
                        currVal--;
                    }
                    if (currVal < minVal) {
                        minVal = currVal;
                        indexMinVal = l + c;
                    }
                    precomputedRMQSmallBlocks[mask][l].push_back(indexMinVal); //[mask][l][l+c]
                }
            }
        }
    }

    void LCA::populateBlockArrays() {
        for (index l = 0; l < eulerTour.size(); l += blockSize) {
            auto startOfBlock = depth.begin() + l;
            auto endOfBlock = l + blockSize <= eulerTour.size() ? startOfBlock + blockSize : depth.end();
            auto minInBlock = std::min_element(startOfBlock, endOfBlock);
            depthMinOfBlock.push_back(*minInBlock);
            eulerTourIndexRealizingDepthMinOfBlock.push_back(std::distance(depth.begin(), minInBlock));
            // if block is smaller than blockSize, then mask will still work - because we won't use the further elements.
            depthBlockToMask.push_back(getMask(startOfBlock, endOfBlock));
        }
    }

    void LCA::populateSparseTable() {
        sparseTableDepthMinOfBlockRMQ.resize(depthMinOfBlock.size());
        for (index l = 0; l < depthMinOfBlock.size(); l++) sparseTableDepthMinOfBlockRMQ[l].push_back(l);
        for (int size = 2; size < depthMinOfBlock.size(); size *= 2) {
            for (index l = 0; l + size < depthMinOfBlock.size(); l++) {
                index leftHalfArgMin = sparseTableDepthMinOfBlockRMQ[l].back();
                index rightHalfArgMin = sparseTableDepthMinOfBlockRMQ[l + (size / 2)].back();
                index argMinForCurrentSize = depthMinOfBlock[leftHalfArgMin] <= depthMinOfBlock[rightHalfArgMin] ?
                                             leftHalfArgMin : rightHalfArgMin;
                sparseTableDepthMinOfBlockRMQ[l].push_back(argMinForCurrentSize);
            }
        }
    }

    LCA::LCA(const FullBranchingTree& fbt) : fbt(fbt) {
        representative.assign(fbt.getTree().upperNodeIdBound(), NetworKit::none);
        eulerTourDFS(fbt.getRoot(), 1);
#ifdef MST_DEBUG
        /*
         * From the paper: "The Euler Tour of T is the sequence of nodes we obtain if we write down the label of each
         * node each time it is visited during a DFS. The array of the Euler tour has length 2n-1 because we start at
         * the root and subsequently output a node each time we traverse an edge.
         * We traverse each of the n-1 edges twice, once in each direction."
         */
        assert(eulerTour.size() == 2 * fbt.getTree().numberOfNodes() - 1);
#endif
        blockSize = MST::log2(eulerTour.size()) / 2;
        if (blockSize >= MIN_ALLOWED_BLOCK_SIZE) {
            populatePrecomputedRMQSmallBlocks();
            populateBlockArrays();
            populateSparseTable();
        }
    }


    node LCA::lcaNaive(node u, node v) const {
        std::unordered_set<node> visited;
        std::optional<node> vertex = {u};
        while (vertex.has_value()) {
            visited.insert(vertex.value());
            vertex = fbt.getParent(vertex.value());
        }
        vertex = {v};
        while (!visited.contains(vertex.value())) {
            vertex = fbt.getParent(vertex.value());
        }
        return vertex.value();
    }

    node LCA::lca(node u, node v) const {
        if (blockSize < MIN_ALLOWED_BLOCK_SIZE) {
            return lcaNaive(u, v);
        }
        const auto blocksInfos = QueryToBlockInfos(this, u, v);
        index rmqDepthLR = [&]() {
            if (blocksInfos.theSameBlock()) {
                return querySmallBlock(blocksInfos.left.blockNumber, blocksInfos.left.indexWithinBlock,
                                       blocksInfos.right.indexWithinBlock);
            } else {
                index index1 = querySmallBlock(
                        blocksInfos.left.blockNumber, blocksInfos.left.indexWithinBlock, blockSize - 1);
                index index2 = querySmallBlock(blocksInfos.right.blockNumber, 0, blocksInfos.right.indexWithinBlock);
                index argMin = depth[index1] <= depth[index2] ? index1 : index2;
                if (blocksInfos.adjacentBlocks()) {
                    return argMin;
                } else {
                    index blockL = blocksInfos.left.blockNumber + 1;
                    index blockR = blocksInfos.right.blockNumber - 1;
                    auto k = MST::log2(blockR - blockL);
                    auto twoToK = 1 << k;
                    auto index1ST = sparseTableDepthMinOfBlockRMQ[blockL][k];
                    auto index2ST = sparseTableDepthMinOfBlockRMQ[blockR - twoToK + 1][k];
                    index argMinST = depthMinOfBlock[index1ST] <= depthMinOfBlock[index2ST] ? index1ST : index2ST;
                    if (depth[argMin] <= depthMinOfBlock[argMinST]) {
                        return argMin;
                    } else {
                        return eulerTourIndexRealizingDepthMinOfBlock[argMinST];
                    };
                }
            }
        }();
        return eulerTour[rmqDepthLR];
    }
}
