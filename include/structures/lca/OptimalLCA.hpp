/*
 * OptimalLCA.hpp
 *
 *  Created on: 05.04.2023
 *      Author: Kamil Kropiewnicki
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

inline constexpr NetworKit::count log2(NetworKit::count n) {
    return n > 1 ? 1 + log2(n >> 1) : 0;
}

struct BlockInfo { NetworKit::index block, offset; };

namespace Koala {

/**
 * @ingroup lca
 * <O(n), O(1)> algorithm from Bender & Farach-Colton "The LCA Problem Revisited" (2000).
 */
template <class T>
class OptimalLCA {
 public:
    using depth_t = NetworKit::count;
    using mask_t = NetworKit::index;

    // preprocessing in O(|V|) time
    explicit OptimalLCA(T& tree);

    // query in O(1) time
    [[nodiscard]] NetworKit::node query(NetworKit::node u, NetworKit::node v) const;

    void verify() const;

 protected:
    static const NetworKit::count MIN_ALLOWED_BLOCK_SIZE = 2;
    const T& tree;

    NetworKit::count block_size;  // log(E.size()) / 2
    // E[] in the paper. E[i] is the label of the ith node visited in the Euler tour of T.
    std::vector<NetworKit::node> E;
    // L[] in the paper. Distance from root. L[i] is the depth of node E[i].
    std::vector<depth_t> L;
    // R in the paper; "time-in" in DFS nomenclature. E[R[i]] = i;
    // note that E[] could have many entries for `i` - R[i] points to the first one.
    std::vector<NetworKit::index> R;
    // small_blocks[mask][l][c]: for a given `mask` (+-1 sequence encoded in binary)
    // it stores argMin of <sequence that generated given mask> in interval [l, l+c] (inclusive).
    // small_blocks[mask][0][block_size-1] will give an answer for the whole block.
    std::vector<std::vector<std::vector<NetworKit::index>>> small_blocks;

    // --- BLOCK ARRAYS --- //
    // A' in the paper - A_prim[i] holds minimum element of L[i*block_size...i*(block_size+1)-1]
    std::vector<depth_t> A_prim;
    // Similar to B in the paper. In our implementation, it holds values incremented by offset
    std::vector<NetworKit::index> B;
    // depthBlock may be described as (startingValue, mask), where mask is a +-1s sequence.
    // The vector below provides mapping: L[i*block_size...i*(block_size+1)-1] -> mask
    std::vector<mask_t> depthBlockToMask;

    // --- SPARSE TABLE --- //
    // M in the paper. M[l][j] is the index of minimum element in A_prim[i][l...l+2^j)
    std::vector<std::vector<NetworKit::index>> M;

    void euler_DFS(NetworKit::node u, NetworKit::count level);
    void populate_precomputed_RMQ_small_blocks();
    void populate_block_arrays();
    void populate_sparse_table();

    [[nodiscard]]
    static mask_t get_mask(
        std::vector<depth_t>::iterator begin, std::vector<depth_t>::iterator end);
    [[nodiscard]]
    NetworKit::index query_small_block(mask_t mask, NetworKit::index l, NetworKit::index r) const;
    [[nodiscard]]
    NetworKit::node query_naive(NetworKit::node u, NetworKit::node v) const;
};

template <class T>
OptimalLCA<T>::OptimalLCA(T& tree) : tree(tree) {
    R.assign(tree.getTree().upperNodeIdBound(), NetworKit::none);
    euler_DFS(tree.getRoot(), 1);
    block_size = log2(E.size()) / 2;
    if (block_size >= MIN_ALLOWED_BLOCK_SIZE) {
        populate_precomputed_RMQ_small_blocks();
        populate_block_arrays();
        populate_sparse_table();
    }
}

template <class T>
NetworKit::node OptimalLCA<T>::query(NetworKit::node u, NetworKit::node v) const {
    if (block_size < MIN_ALLOWED_BLOCK_SIZE) {
        return query_naive(u, v);
    }
    auto [l, r] = std::minmax(R[u], R[v]);
    BlockInfo left{l / block_size, l % block_size}, right{r / block_size, r % block_size};
    if (left.block == right.block) {
        return E[query_small_block(left.block, left.offset, right.offset)];
    }
    NetworKit::index index_left = query_small_block(left.block, left.offset, block_size - 1);
    NetworKit::index index_right = query_small_block(right.block, 0, right.offset);
    auto arg_min = L[index_left] <= L[index_right] ? index_left : index_right;
    if (left.block + 1 == right.block) {
        return E[arg_min];
    }
    auto k = log2(right.block - left.block - 2);
    index_left = M[left.block + 1][k], index_right = M[right.block - (1 << k)][k];
    auto arg_min_st = A_prim[index_left] <= A_prim[index_right] ? index_left : index_right;
    return L[arg_min] <= A_prim[arg_min_st] ? E[arg_min] : E[B[arg_min_st]];
}

template <class T>
void OptimalLCA<T>::verify() const {
    // From the paper: "The Euler Tour of T is the sequence of nodes we obtain if we write down
    // the label of each node each time it is visited during a DFS. The array of the Euler tour
    // has length 2n-1 because we start at the root and subsequently output a node each time
    // we traverse an edge. We traverse each of the n-1 edges twice, once in each direction."
    assert(E.size() == 2 * tree.getTree().numberOfNodes() - 1);
}

template <class T>
void OptimalLCA<T>::euler_DFS(const NetworKit::node u, const depth_t level) {
    // We cannot use NetworKit's DFS because it accepts just one handle that is called
    // during traversing adjacency list before following an edge - but we need the second handle
    // to call after coming back from the recursion.
    if (R[u] == NetworKit::none) {
        R[u] = E.size();
    }
    E.push_back(u), L.push_back(level);
    tree.getTree().forEdgesOf(u, [this, level](NetworKit::node u, NetworKit::node child) {
        euler_DFS(child, level + 1);
        E.push_back(u), L.push_back(level);
    });
}

template <class T>
NetworKit::index OptimalLCA<T>::query_small_block(
        NetworKit::index block, NetworKit::index l, NetworKit::index r) const {
    auto mask = depthBlockToMask[block];
    return small_blocks[mask][l][r - l] + (block * block_size);
}

template <class T>
void OptimalLCA<T>::populate_precomputed_RMQ_small_blocks() {
    // `block_size-1` because we count +-1s between elements.
    const NetworKit::count masks_count = 1 << (block_size - 1);
    small_blocks.resize(masks_count, std::vector<std::vector<NetworKit::index>>(block_size));
    for (NetworKit::count mask = 0; mask < masks_count; mask++) {  // O(sqrt(|V|)) iterations
        for (NetworKit::index l = 0; l < block_size; l++) {  // O(log(|V|)) iterations
            small_blocks[mask][l].push_back(l);
            int minimum = 0, current = 0;
            NetworKit::index minimum_index = l;
            // c - how many +-1 are considered. corresponds to mask[l...l+c]
            for (size_t c = 1; l + c < block_size; c++) {  // O(log(|V|))) iterations
                // skip first `l` elements and include `c` +-1s looking to the right.
                if (mask & (1 << (l + c - 1))) {
                    current++;
                } else {
                    current--;
                }
                if (current < minimum) {
                    minimum = current;
                    minimum_index = l + c;
                }
                small_blocks[mask][l].push_back(minimum_index);  // set [mask][l][l+c]
            }
        }
    }
}

template <class T>
void OptimalLCA<T>::populate_block_arrays() {
    for (NetworKit::index l = 0; l < E.size(); l += block_size) {
        auto startOfBlock = L.begin() + l;
        auto endOfBlock = l + block_size <= E.size() ? startOfBlock + block_size : L.end();
        auto minInBlock = std::min_element(startOfBlock, endOfBlock);
        A_prim.push_back(*minInBlock);
        B.push_back(std::distance(L.begin(), minInBlock));
        // if block is smaller than block_size, then mask will still work ignoring further elements
        depthBlockToMask.push_back(get_mask(startOfBlock, endOfBlock));
    }
}

// Mask: traverse input from the left, but construct int representing mask from the right.
// In mask, LeastSignificantBit tells us about the first difference in a sequence.
template <class T>
OptimalLCA<T>::mask_t OptimalLCA<T>::get_mask(
        std::vector<OptimalLCA<T>::depth_t>::iterator begin,
        std::vector<OptimalLCA<T>::depth_t>::iterator end) {
    mask_t mask = 0;
    std::advance(begin, 1);
    for (NetworKit::index i = 0; begin != end; i++, std::advance(begin, 1)) {
        mask |= ((*begin) > (*prev(begin))) << i;
    }
    return mask;
}

template <class T>
void OptimalLCA<T>::populate_sparse_table() {
    M.resize(A_prim.size());
    for (NetworKit::index l = 0; l < A_prim.size(); l++) {
        M[l].push_back(l);
    }
    for (int size = 2; size < A_prim.size(); size *= 2) {
        for (NetworKit::index l = 0; l + size < A_prim.size(); l++) {
            NetworKit::index leftArgMin = M[l].back(), rightArgMin = M[l + (size / 2)].back();
            NetworKit::index argMinForCurrentSize =
                A_prim[leftArgMin] <= A_prim[rightArgMin] ? leftArgMin : rightArgMin;
            M[l].push_back(argMinForCurrentSize);
        }
    }
}

template <class T>
NetworKit::node OptimalLCA<T>::query_naive(NetworKit::node u, NetworKit::node v) const {
    std::unordered_set<NetworKit::node> visited;
    std::optional<NetworKit::node> vertex = {u};
    while (vertex.has_value()) {
        visited.insert(vertex.value());
        vertex = tree.getParent(vertex.value());
    }
    vertex = {v};
    while (!visited.contains(vertex.value())) {
        vertex = tree.getParent(vertex.value());
    }
    return vertex.value();
}

}  // namespace Koala
