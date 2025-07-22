#include "shortest_path/planar/TopologyBasedHeap.hpp"

#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <vector>

namespace Koala {
static const NetworKit::edgeweight INF = std::numeric_limits<NetworKit::edgeweight>::max();

// root is stored in storage[1] to simplify traverse logic
void TopologyHeap::fix_index(NetworKit::index i) {
    int l = i * 2, r = i * 2 + 1;

    if (storage[l].first < storage[r].first) {
        storage[i] = storage[l];
    } else {
        storage[i] = storage[r];
    }
}

void TopologyHeap::initial_distances() {
    // initial distances to boundary nodes from source

    if (node_regions[source].size() > 1 ||
        (extra_boundary_nodes.find(source) !=
            extra_boundary_nodes.end())) {  // boundary node distances are already calculated
        storage[node_storage_index[source]].first = 0;
    }

    auto initial_region_id = node_regions[source][0];

    auto& region = regions[initial_region_id];
    auto subgraph = NetworKit::GraphTools::subgraphFromNodes(
        graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
    std::vector<NetworKit::node> boundary;

    for (auto batch_id : region_batches[initial_region_id]) {
        for (auto b : batches[batch_id]) {
            boundary.push_back(b);
        }
    }
    NetworKit::MultiTargetDijkstra dijkstra(subgraph, source, boundary.begin(), boundary.end());
    dijkstra.run();
    auto map = dijkstra.getTargetIndexMap();
    auto distances = dijkstra.getDistances();
    for (auto b : boundary) {
        storage[node_storage_index[b]].first = distances[map[b]];
    }
}

void TopologyHeap::initialize_storage() {
    // initialize the heap to have boundary sets in clusters
    std::vector<NetworKit::index> last_batch;
    region_batches.resize(regions.size());
    for (NetworKit::index i = 0; i < sorted_boundary_node_regions.size(); i++) {
        NetworKit::index node_index = sorted_boundary_node_regions[i].back();
        sorted_boundary_node_regions[i].pop_back();

        storage[i + storage_size] = {INF, node_index};
        node_storage_index[node_index] = i + storage_size;

        if (last_batch == sorted_boundary_node_regions[i]) {
            batches.back().push_back(node_index);
        } else {
            last_batch = sorted_boundary_node_regions[i];
            batches.push_back({node_index});
            for (auto region : sorted_boundary_node_regions[i]) {
                region_batches[region].push_back(batches.size() - 1);
            }
        }
    }

    initial_distances();

    // build the heap
    for (NetworKit::index i = storage_size - 1; i >= 1; i--) {
        fix_index(i);
    }
}

TopologyHeap::TopologyHeap(NetworKit::Graph& graph, node_subsets_t& regions,
    pair_distance_t& distances, NetworKit::node source,
    const std::unordered_set<NetworKit::node>& extra_boundary_nodes)
    : graph(graph),
      regions(regions),
      distances(distances),
      source(source),
      extra_boundary_nodes(extra_boundary_nodes) {
    // get list per node of regions a node is a part of
    node_regions.assign(graph.numberOfNodes(), std::vector<NetworKit::index>{});
    node_storage_index.resize(graph.numberOfNodes());
    is_closed.assign(graph.numberOfNodes(), 0);
    for (NetworKit::index i = 0; i < regions.size(); i++) {
        for (auto& node : regions[i]) {
            node_regions[node].push_back(i);
        }
    }

    // filter boundary nodes and add node_id to the end
    for (NetworKit::index i = 0; i < graph.numberOfNodes(); i++) {
        if (node_regions[i].size() > 1) {
            node_regions[i].push_back(i);
            sorted_boundary_node_regions.push_back(node_regions[i]);
            node_regions[i].pop_back();
        } else if (extra_boundary_nodes.find(i) != extra_boundary_nodes.end()) {
            sorted_boundary_node_regions.push_back({node_regions[i][0], i});
        }
    }

    size = sorted_boundary_node_regions.size();
    std::sort(sorted_boundary_node_regions.begin(), sorted_boundary_node_regions.end());
    // TODO(289Adam289): change sort from std::sort to some linear equivalent

    storage_size = 1;
    while (storage_size < size) storage_size <<= 1;

    // double the storage_size. Actual storage and heap structure
    storage.assign(storage_size << 1, {INF, NetworKit::none});

    initialize_storage();
}

std::pair<NetworKit::edgeweight, NetworKit::node> TopologyHeap::top() {
    if (size > 0) {
        assert(storage[1].second != NetworKit::none);
    }
    return storage[1];
}

void TopologyHeap::batch_update(
    NetworKit::index batch_id, NetworKit::node node, NetworKit::edgeweight current_value) {
    auto batch = batches[batch_id];
    NetworKit::index node_id = node_storage_index[node];
    std::deque<NetworKit::index> index_to_fix;
    for (auto b : batch) {
        auto storage_id = node_storage_index[b];
        if (index_to_fix.empty() || index_to_fix.back() != storage_id >> 1) {
            index_to_fix.emplace_back(storage_id >> 1);
        }
        if (is_closed[b]) continue;
        if (b == node) {
            continue;
        } else {
            auto distance = INF;

            if (distances.find({node, b}) != distances.end()) {
                distance = current_value + distances[{node, b}];
            }

            storage[storage_id].first = std::min(storage[storage_id].first, distance);
        }
    }

    while (!index_to_fix.empty()) {
        auto id = index_to_fix.front();
        index_to_fix.pop_front();

        if (id == 0) break;

        fix_index(id);

        if (index_to_fix.back() != id >> 1) {
            index_to_fix.emplace_back(id >> 1);
        }
    }
}

void TopologyHeap::close(NetworKit::node node) {
    size--;
    assert(node_regions[node].size() > 1 ||
           (extra_boundary_nodes.find(source) != extra_boundary_nodes.end()));

    auto current_value = storage[node_storage_index[node]].first;
    storage[node_storage_index[node]].first = INF;
    is_closed[node] = true;
    std::unordered_set<NetworKit::index> batches_to_update;
    for (auto regionId : node_regions[node]) {
        for (auto batch_id : region_batches[regionId]) {
            batches_to_update.insert(batch_id);
        }
    }
    for (auto id : batches_to_update) {
        batch_update(id, node, current_value);
    }
}

bool TopologyHeap::empty() {
    return size <= 0;
}

}  // namespace Koala
