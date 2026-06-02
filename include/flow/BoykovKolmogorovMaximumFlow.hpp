#pragma once

#include <queue>
#include <unordered_map>

#include "flow/MaximumFlow.hpp"

namespace Koala {

/**
 * @ingroup flow
 * Boykov-Kolmogorov augmenting-path maximum-flow algorithm.
 *
 * The algorithm grows source and sink search trees, augments when the trees meet, and adopts
 * orphaned vertices after saturated tree edges are removed.
 *
 * @see https://doi.org/10.1109/TPAMI.2004.60
 */
class BoykovKolmogorovMaximumFlow final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;

    /**
     * Compute a maximum flow with the grow, augment, and adopt phases.
     */
    void run();

 private:
    enum class NodeType { SOURCE, TARGET, FREE };

    std::unordered_map<NetworKit::Edge, int> flow;
    NetworKit::node spath, tpath;
    std::unordered_map<NetworKit::node, NetworKit::node> parent;
    std::unordered_map<NetworKit::node, NodeType> tree;
    std::unordered_map<NetworKit::Edge, int> capacity;
    std::queue<NetworKit::node> active, orphan;

    int tree_capacity(NetworKit::node, NetworKit::node);
    void initialize();
    bool grow();
    int augment();
    void adopt();
    bool origin(NetworKit::node);
};

}  // namespace Koala
