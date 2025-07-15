#include <set>

#include "shortest_path/HenzingerPlanarSSSP.hpp"

namespace Koala {

static const int INF = std::numeric_limits<int>::max();

// Updates key of item under id. In case set does not contain the item inserts it.
void HenzingerPlanarSSSP::HenzingerPriorityQueue::update(
    NetworKit::node id, NetworKit::count new_key) {
    if (id_map.find(id) == id_map.end()) {
        push(id, new_key);
        return;
    }

    int old_key = id_map[id];
    assert(set.erase({old_key, id}) == 1 && "update");

    id_map[id] = new_key;
    set.insert({new_key, id});

    minimum_element = *(set.begin());
}
void HenzingerPlanarSSSP::HenzingerPriorityQueue::push(NetworKit::node id, NetworKit::count key) {
    assert(id_map.find(id) == id_map.end() && "insert");

    set.insert({key, id});
    id_map[id] = key;

    minimum_element = *(set.begin());
}
void HenzingerPlanarSSSP::HenzingerPriorityQueue::deactivate(NetworKit::node id) {
    assert(id_map.find(id) != id_map.end() && "deactivate");

    int old_key = id_map[id];
    assert(set.erase({old_key, id}) == 1);
    id_map.erase(id);

    if (!empty()) {
        minimum_element = *(set.begin());
    }
}
NetworKit::node HenzingerPlanarSSSP::HenzingerPriorityQueue::minimum_item() {
    assert(!empty() && "minimum_item");
    return minimum_element.second;
}
NetworKit::count HenzingerPlanarSSSP::HenzingerPriorityQueue::minimum_key() {
    assert(!empty() && "minimum_key");
    return minimum_element.first;
}
bool HenzingerPlanarSSSP::HenzingerPriorityQueue::empty() {
    return set.size() == 0;
}
}  // namespace Koala
