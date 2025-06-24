#pragma once

#include <networkit/graph/Graph.hpp>
#include <set>

static const int INF = std::numeric_limits<int>::max();

namespace Koala {
class HPriorityQueue {
    std::set<std::pair<NetworKit::count, NetworKit::node>> set;
    std::unordered_map<NetworKit::node, NetworKit::count> idMap;

    // caching the mininmal elemnet to ensure constan time read
    std::pair<NetworKit::count, NetworKit::node> minElement;

   public:
    void initialize(std::vector<NetworKit::node>& ids);
    void insert(NetworKit::node id, NetworKit::count key);
    void updateKey(NetworKit::node id, NetworKit::count newKey);
    void deactivateItem(NetworKit::node id);
    NetworKit::count minKey();
    NetworKit::node minItem();
    bool empty();
};

} /* namespace Koala */
