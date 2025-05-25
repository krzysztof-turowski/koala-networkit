#include "planar_sssp/HenzingerPriorityQueue.hpp"


namespace Koala {
    //requires for ids to be unique
    void HPriorityQueue::initialize(std::vector<NetworKit::node>& ids) {
        for (auto id : ids) {
            set.insert({ INF, id });
            idMap[id] = INF;
        }

        if (!empty()) {
            minElement = *(set.begin());
        }
    }

    //Updates key of item under id. In case set does not contain the item inserts it.
    void HPriorityQueue::updateKey(NetworKit::node id, NetworKit::count newKey) {
        if (idMap.find(id) == idMap.end()) {
            insert(id, newKey);
            return;
        }

        int oldKey = idMap[id];
        assert(set.erase({ oldKey, id }) == 1 && "updateKey");

        idMap[id] = newKey;
        set.insert({ newKey, id });

        minElement = *(set.begin());
    }
    void HPriorityQueue::insert(NetworKit::node id, NetworKit::count key) {
        assert(idMap.find(id) == idMap.end() && "insert");

        set.insert({ key, id });
        idMap[id] = key;

        minElement = *(set.begin());
    }
    void HPriorityQueue::deactivateItem(NetworKit::node id) {
        assert(idMap.find(id) != idMap.end() && "deactivateItem");

        int oldKey = idMap[id];
        assert(set.erase({ oldKey, id }) == 1);
        idMap.erase(id);

        if (!empty()) {
            minElement = *(set.begin());
        }
    }
    NetworKit::node HPriorityQueue::minItem() {
        assert(!empty() && "minItem");
        return minElement.second;
    }
    NetworKit::count HPriorityQueue::minKey() {
        assert(!empty() && "minKey");
        return minElement.first;
    }
    bool HPriorityQueue::empty() { return set.size() == 0; }
}