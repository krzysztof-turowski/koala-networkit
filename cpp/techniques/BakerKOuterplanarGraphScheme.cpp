#include <algorithm>
#include <cstddef>
#include <deque>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <networkit/components/ConnectedComponents.hpp>

#include "techniques/BakerKOuterplanarGraphScheme.hpp"

namespace {

constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();

struct IndexPairHash {
    std::size_t operator()(
            const std::pair<std::size_t, std::size_t> &value) const noexcept {
        const auto first = std::hash<std::size_t>{}(value.first);
        const auto second = std::hash<std::size_t>{}(value.second);
        return first ^ (second << 1);
    }
};

class DisjointSet {
 public:
    explicit DisjointSet(std::size_t) { }

    std::size_t find(std::size_t element) {
        const auto [iterator, inserted] =
            parent_.try_emplace(element, element);
        rank_.try_emplace(element, 0);
        if (!inserted && iterator->second != element) {
            iterator->second = find(iterator->second);
        }
        return iterator->second;
    }

    bool unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first == second) {
            return false;
        }
        if (rank_[first] < rank_[second]) {
            std::swap(first, second);
        }
        parent_[second] = first;
        if (rank_[first] == rank_[second]) {
            ++rank_[first];
        }
        return true;
    }

 private:
    std::unordered_map<std::size_t, std::size_t> parent_;
    std::unordered_map<std::size_t, unsigned char> rank_;
};

class PlaneGraph {
 public:
    struct Dart {
        NetworKit::node from;
        NetworKit::node to;
        std::size_t twin;
        std::size_t edge;
        std::size_t rotationPrevious = invalid_index;
        std::size_t rotationNext = invalid_index;
    };

    struct Edge {
        std::size_t first;
        std::size_t second;
        bool original;
        std::optional<Koala::BakerFakeEdgeType> fakeType;
    };

    struct Faces {
        std::vector<std::vector<std::size_t>> boundaries;
        std::vector<std::size_t> faceOfDart;
    };

    explicit PlaneGraph(std::size_t node_bound) noexcept
        : nodeBound_(node_bound) { }

    std::size_t addDetachedEdge(
            NetworKit::node u, NetworKit::node v, bool original,
            std::optional<Koala::BakerFakeEdgeType> fake_type) {
        const std::size_t edge_index = edges_.size();
        const std::size_t first = darts_.size();
        const std::size_t second = first + 1;
        darts_.push_back({u, v, second, edge_index});
        darts_.push_back({v, u, first, edge_index});
        edges_.push_back({first, second, original, fake_type});
        adjacency_[u].insert(v);
        adjacency_[v].insert(u);
        return edge_index;
    }

    void setRotation(
            NetworKit::node vertex,
            const std::vector<std::size_t> &rotation) {
        if (rotation.empty()) {
            anchor_.erase(vertex);
            return;
        }
        anchor_[vertex] = rotation.front();
        for (std::size_t i = 0; i < rotation.size(); ++i) {
            const std::size_t current = rotation[i];
            if (darts_[current].from != vertex) {
                throw std::logic_error("A plane rotation contains a foreign dart");
            }
            darts_[current].rotationPrevious =
                rotation[(i + rotation.size() - 1) % rotation.size()];
            darts_[current].rotationNext =
                rotation[(i + 1) % rotation.size()];
        }
    }

    std::pair<std::size_t, std::size_t> addEdgeInFace(
            std::size_t first_incoming, std::size_t first_outgoing,
            std::size_t second_incoming, std::size_t second_outgoing,
            bool original,
            std::optional<Koala::BakerFakeEdgeType> fake_type) {
        const auto first_vertex = darts_[first_outgoing].from;
        const auto second_vertex = darts_[second_outgoing].from;
        if (first_vertex == second_vertex) {
            throw std::logic_error("A Baker triangulation attempted a loop");
        }
        if (darts_[first_incoming].to != first_vertex
            || darts_[second_incoming].to != second_vertex
            || faceNext(first_incoming) != first_outgoing
            || faceNext(second_incoming) != second_outgoing) {
            throw std::logic_error("Invalid face wedges for a Baker edge");
        }

        const std::size_t edge_index = addDetachedEdge(
            first_vertex, second_vertex, original, fake_type);
        const std::size_t first_dart = edges_[edge_index].first;
        const std::size_t second_dart = edges_[edge_index].second;
        insertInFaceWedge(
            first_incoming, first_outgoing, first_dart);
        insertInFaceWedge(
            second_incoming, second_outgoing, second_dart);
        return {first_dart, second_dart};
    }

    std::size_t faceNext(std::size_t dart) const {
        return darts_[darts_[dart].twin].rotationPrevious;
    }

    Faces faces() const {
        Faces result;
        result.faceOfDart.assign(darts_.size(), invalid_index);
        for (std::size_t start = 0; start < darts_.size(); ++start) {
            if (result.faceOfDart[start] != invalid_index) {
                continue;
            }
            const std::size_t face = result.boundaries.size();
            std::vector<std::size_t> boundary;
            std::size_t dart = start;
            do {
                if (result.faceOfDart[dart] != invalid_index) {
                    throw std::logic_error("A plane face traversal is inconsistent");
                }
                result.faceOfDart[dart] = face;
                boundary.push_back(dart);
                dart = faceNext(dart);
            } while (dart != start);
            result.boundaries.push_back(std::move(boundary));
        }
        return result;
    }

    std::vector<std::size_t> rotation(NetworKit::node vertex) const {
        std::vector<std::size_t> result;
        const auto anchor = anchor_.find(vertex);
        if (anchor == anchor_.end()) {
            return result;
        }
        const std::size_t start = anchor->second;
        std::size_t dart = start;
        do {
            result.push_back(dart);
            dart = darts_[dart].rotationNext;
        } while (dart != start);
        return result;
    }

    bool hasEdge(NetworKit::node u, NetworKit::node v) const {
        const auto neighbors = adjacency_.find(u);
        return neighbors != adjacency_.end()
            && neighbors->second.contains(v);
    }

    const Dart &dart(std::size_t index) const {
        return darts_[index];
    }

    Dart &dart(std::size_t index) {
        return darts_[index];
    }

    const Edge &edge(std::size_t index) const {
        return edges_[index];
    }

    const std::vector<Edge> &edges() const {
        return edges_;
    }

    std::size_t numberOfDarts() const {
        return darts_.size();
    }

    std::size_t nodeBound() const {
        return nodeBound_;
    }

 private:
    void insertInFaceWedge(
            std::size_t incoming, std::size_t outgoing,
            std::size_t inserted) {
        const std::size_t before = darts_[incoming].twin;
        if (darts_[before].rotationPrevious != outgoing) {
            throw std::logic_error("A Baker edge was inserted outside its face");
        }
        if (before == outgoing) {
            darts_[before].rotationPrevious = inserted;
            darts_[before].rotationNext = inserted;
            darts_[inserted].rotationPrevious = before;
            darts_[inserted].rotationNext = before;
            return;
        }
        darts_[before].rotationPrevious = inserted;
        darts_[inserted].rotationNext = before;
        darts_[inserted].rotationPrevious = outgoing;
        darts_[outgoing].rotationNext = inserted;
    }

    std::vector<Dart> darts_;
    std::vector<Edge> edges_;
    std::unordered_map<NetworKit::node, std::size_t> anchor_;
    std::unordered_map<
        NetworKit::node, std::unordered_set<NetworKit::node>> adjacency_;
    std::size_t nodeBound_;
};

using BoostGraph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
    boost::property<boost::edge_index_t, std::size_t>>;
using BoostEdge = boost::graph_traits<BoostGraph>::edge_descriptor;
using Embedding = std::vector<std::vector<BoostEdge>>;
using LevelMap =
    std::unordered_map<NetworKit::node, NetworKit::count>;

void assignEdgeIndices(BoostGraph &graph) {
    auto edge_index = get(boost::edge_index, graph);
    std::size_t index = 0;
    auto [iterator, end] = edges(graph);
    for (; iterator != end; ++iterator) {
        put(edge_index, *iterator, index);
        ++index;
    }
}

BoostGraph makeBoostGraph(
        const NetworKit::Graph &graph,
        const std::vector<NetworKit::node> &component,
        std::vector<NetworKit::node> &local_to_global) {
    local_to_global = component;
    std::unordered_map<NetworKit::node, std::size_t> global_to_local;
    for (std::size_t i = 0; i < local_to_global.size(); ++i) {
        global_to_local[local_to_global[i]] = i;
    }

    BoostGraph result(local_to_global.size());
    for (std::size_t local_u = 0;
         local_u < local_to_global.size(); ++local_u) {
        const auto global_u = local_to_global[local_u];
        graph.forNeighborsOf(
            global_u, [&](NetworKit::node global_v) {
                if (global_u >= global_v) {
                    return;
                }
                const auto local_v = global_to_local.find(global_v);
                if (local_v != global_to_local.end()) {
                    add_edge(local_u, local_v->second, result);
                }
            });
    }
    assignEdgeIndices(result);
    return result;
}

PlaneGraph makePlaneGraph(
        const BoostGraph &boost_graph, const Embedding &embedding,
        const std::vector<NetworKit::node> &local_to_global,
        std::size_t node_bound) {
    PlaneGraph result(node_bound);
    std::vector<std::size_t> boost_to_plane(num_edges(boost_graph), invalid_index);
    auto edge_index = get(boost::edge_index, boost_graph);
    auto [iterator, end] = edges(boost_graph);
    for (; iterator != end; ++iterator) {
        const auto local_u = source(*iterator, boost_graph);
        const auto local_v = target(*iterator, boost_graph);
        const std::size_t plane_edge = result.addDetachedEdge(
            local_to_global[local_u], local_to_global[local_v], true,
            std::nullopt);
        boost_to_plane[get(edge_index, *iterator)] = plane_edge;
    }

    for (std::size_t local = 0; local < embedding.size(); ++local) {
        const auto global = local_to_global[local];
        std::vector<std::size_t> rotation;
        rotation.reserve(embedding[local].size());
        for (const auto boost_edge : embedding[local]) {
            const auto &plane_edge =
                result.edge(boost_to_plane[get(edge_index, boost_edge)]);
            const auto first = result.dart(plane_edge.first);
            rotation.push_back(
                first.from == global ? plane_edge.first : plane_edge.second);
        }
        result.setRotation(global, rotation);
    }
    return result;
}

Embedding restrictEmbedding(
        const BoostGraph &boost_graph,
        const std::vector<NetworKit::node> &local_to_global,
        const Koala::BakerPlaneEmbedding &embedding) {
    std::unordered_map<NetworKit::node, std::size_t> global_to_local;
    global_to_local.reserve(local_to_global.size());
    for (std::size_t local = 0; local < local_to_global.size(); ++local) {
        global_to_local.emplace(local_to_global[local], local);
    }

    Embedding result(local_to_global.size());
    for (std::size_t local = 0; local < local_to_global.size(); ++local) {
        const auto global = local_to_global[local];
        if (global >= embedding.rotation.size()) {
            throw std::invalid_argument(
                "A prescribed Baker embedding has no vertex rotation");
        }
        std::unordered_set<std::size_t> seen_neighbors;
        for (const auto global_neighbor : embedding.rotation[global]) {
            const auto local_neighbor = global_to_local.find(global_neighbor);
            if (local_neighbor == global_to_local.end()) {
                throw std::invalid_argument(
                    "A prescribed Baker rotation contains a foreign vertex");
            }
            const auto [boost_edge, exists] = edge(
                local, local_neighbor->second, boost_graph);
            if (!exists || !seen_neighbors.insert(local_neighbor->second).second) {
                throw std::invalid_argument(
                    "A prescribed Baker rotation contains an invalid edge");
            }
            result[local].push_back(boost_edge);
        }
        if (result[local].size() != degree(local, boost_graph)) {
            throw std::invalid_argument(
                "A prescribed Baker rotation omits a graph edge");
        }
    }
    return result;
}

std::size_t selectOuterFace(
        const PlaneGraph &, const PlaneGraph::Faces &faces) {
    std::size_t best_face = 0;
    std::size_t best_length = 0;
    for (std::size_t face = 0; face < faces.boundaries.size(); ++face) {
        if (faces.boundaries[face].size() > best_length) {
            best_face = face;
            best_length = faces.boundaries[face].size();
        }
    }
    return best_face;
}

std::size_t selectOuterFace(
        const PlaneGraph &graph, const PlaneGraph::Faces &faces,
        std::vector<NetworKit::node> requested_vertices) {
    std::sort(requested_vertices.begin(), requested_vertices.end());
    requested_vertices.erase(
        std::unique(
            requested_vertices.begin(), requested_vertices.end()),
        requested_vertices.end());

    for (std::size_t face = 0; face < faces.boundaries.size(); ++face) {
        std::vector<NetworKit::node> boundary_vertices;
        boundary_vertices.reserve(faces.boundaries[face].size());
        for (const auto dart : faces.boundaries[face]) {
            boundary_vertices.push_back(graph.dart(dart).from);
        }
        std::sort(boundary_vertices.begin(), boundary_vertices.end());
        boundary_vertices.erase(
            std::unique(
                boundary_vertices.begin(), boundary_vertices.end()),
            boundary_vertices.end());
        if (boundary_vertices == requested_vertices) {
            return face;
        }
    }
    throw std::invalid_argument(
        "Requested Baker outer-face vertices do not form an embedding face");
}

LevelMap computeLevels(
        const PlaneGraph &graph, const PlaneGraph::Faces &faces,
        std::size_t outer_face) {
    const std::size_t infinity = std::numeric_limits<std::size_t>::max();
    std::unordered_map<
        NetworKit::node, std::vector<std::size_t>> incident_faces;
    for (std::size_t face = 0; face < faces.boundaries.size(); ++face) {
        std::unordered_set<NetworKit::node> seen;
        for (const auto dart : faces.boundaries[face]) {
            const auto vertex = graph.dart(dart).from;
            if (seen.insert(vertex).second) {
                incident_faces[vertex].push_back(face);
            }
        }
    }

    struct QueueEntry {
        bool face;
        std::size_t index;
    };
    std::vector<std::size_t> face_distance(
        faces.boundaries.size(), infinity);
    std::unordered_map<NetworKit::node, std::size_t> vertex_distance;
    std::deque<QueueEntry> queue;
    face_distance[outer_face] = 0;
    queue.push_back({true, outer_face});

    while (!queue.empty()) {
        const auto current = queue.front();
        queue.pop_front();
        if (current.face) {
            const auto next_distance = face_distance[current.index] + 1;
            for (const auto dart : faces.boundaries[current.index]) {
                const auto vertex = graph.dart(dart).from;
                const bool inserted =
                    vertex_distance.emplace(
                        vertex, next_distance).second;
                if (inserted) {
                    queue.push_back({false, vertex});
                }
            }
        } else {
            const auto next_distance =
                vertex_distance.at(current.index) + 1;
            for (const auto face : incident_faces.at(current.index)) {
                if (face_distance[face] == infinity) {
                    face_distance[face] = next_distance;
                    queue.push_back({true, face});
                }
            }
        }
    }

    LevelMap result;
    result.reserve(vertex_distance.size());
    for (const auto &[vertex, distance] : vertex_distance) {
        result.emplace(
            vertex,
            static_cast<NetworKit::count>((distance + 1) / 2));
    }
    return result;
}

std::pair<NetworKit::node, NetworKit::node> orderedEdge(
        NetworKit::node u, NetworKit::node v) {
    return u < v ? std::make_pair(u, v) : std::make_pair(v, u);
}

void triangulateSimpleBoundary(
        PlaneGraph &graph, const std::vector<std::size_t> &boundary,
        const LevelMap &levels,
        DisjointSet &same_level_components,
        std::vector<Koala::BakerFakeEdge> &fake_edges) {
    if (boundary.size() <= 3) {
        return;
    }

    struct BoundaryEntry {
        std::size_t dart;
        std::size_t previous;
        std::size_t next;
        bool active = true;
        bool queued = false;
    };

    std::vector<BoundaryEntry> entries;
    entries.reserve(boundary.size());
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        entries.push_back({
            boundary[i],
            (i + boundary.size() - 1) % boundary.size(),
            (i + 1) % boundary.size()});
    }
    std::size_t active_count = entries.size();

    auto is_eligible = [&](std::size_t index) {
        if (!entries[index].active || active_count <= 3) {
            return false;
        }
        const auto previous = entries[index].previous;
        const auto u = graph.dart(entries[previous].dart).from;
        const auto v = graph.dart(entries[index].dart).to;
        const auto difference = levels.at(u) > levels.at(v)
            ? levels.at(u) - levels.at(v)
            : levels.at(v) - levels.at(u);
        return u != v && difference <= 1
            && (difference != 0
                || same_level_components.find(u)
                    != same_level_components.find(v));
    };

    std::deque<std::size_t> ears;
    auto enqueue = [&](std::size_t index) {
        if (!entries[index].queued && is_eligible(index)) {
            entries[index].queued = true;
            ears.push_back(index);
        }
    };
    for (std::size_t i = 0; i < entries.size(); ++i) {
        enqueue(i);
    }

    while (!ears.empty() && active_count > 3) {
        const auto ear = ears.front();
        ears.pop_front();
        entries[ear].queued = false;
        if (!is_eligible(ear)) {
            continue;
        }

        const auto previous = entries[ear].previous;
        const auto previous_previous = entries[previous].previous;
        const auto next = entries[ear].next;
        const auto u = graph.dart(entries[previous].dart).from;
        const auto v = graph.dart(entries[ear].dart).to;
        Koala::BakerFakeEdgeType type =
            Koala::BakerFakeEdgeType::TRIANGULATION;
        if (levels.at(u) == levels.at(v)
            && same_level_components.unite(u, v)) {
            type = Koala::BakerFakeEdgeType::CONNECTIVITY;
        }

        const auto [forward, backward] = graph.addEdgeInFace(
            entries[previous_previous].dart, entries[previous].dart,
            entries[ear].dart, entries[next].dart, false, type);
        static_cast<void>(backward);
        const auto [edge_u, edge_v] = orderedEdge(u, v);
        fake_edges.push_back({edge_u, edge_v, type});

        entries[previous].dart = forward;
        entries[previous].next = next;
        entries[next].previous = previous;
        entries[ear].active = false;
        --active_count;
        enqueue(previous);
        enqueue(next);
    }
}

void triangulateBoundary(
        PlaneGraph &graph, const std::vector<std::size_t> &boundary,
        const LevelMap &levels,
        DisjointSet &same_level_components,
        std::vector<Koala::BakerFakeEdge> &fake_edges) {
    triangulateSimpleBoundary(
        graph, boundary, levels,
        same_level_components, fake_edges);
}

void triangulateInternalFaces(
        PlaneGraph &graph, std::size_t outer_face,
        const PlaneGraph::Faces &initial_faces,
        const LevelMap &levels,
        std::vector<Koala::BakerFakeEdge> &fake_edges) {
    DisjointSet same_level_components(graph.nodeBound());
    for (const auto &edge : graph.edges()) {
        const auto u = graph.dart(edge.first).from;
        const auto v = graph.dart(edge.first).to;
        if (levels.at(u) == levels.at(v)) {
            same_level_components.unite(u, v);
        }
    }

    for (std::size_t face = 0;
         face < initial_faces.boundaries.size(); ++face) {
        if (face == outer_face) {
            continue;
        }
        triangulateBoundary(
            graph, initial_faces.boundaries[face], levels,
            same_level_components, fake_edges);
    }
}

struct WorkTreeVertex {
    bool face = false;
    std::size_t localFace = invalid_index;
    std::size_t leafDart = invalid_index;
    std::vector<std::size_t> neighbors;
    std::optional<std::size_t> enclosedWork;
};

struct ComponentWork {
    NetworKit::count level = 0;
    std::vector<NetworKit::node> vertices;
    PlaneGraph embedding{0};
    std::vector<std::optional<std::size_t>> globalDart;
    std::vector<std::optional<std::size_t>> referenceGlobalDart;
    std::vector<WorkTreeVertex> treeVertices;
    std::vector<std::size_t> outerLeaves;
    std::vector<std::size_t> outerLeafDarts;
    std::optional<std::size_t> parentWork;
    std::optional<std::size_t> parentFaceVertex;
    std::optional<std::size_t> publicTree;
};

std::vector<std::size_t> findBridges(
        const PlaneGraph &graph,
        const std::vector<NetworKit::node> &vertices) {
    std::unordered_map<NetworKit::node, std::size_t> discovery;
    std::unordered_map<NetworKit::node, std::size_t> low;
    discovery.reserve(vertices.size());
    low.reserve(vertices.size());
    for (const auto vertex : vertices) {
        discovery.emplace(vertex, 0);
        low.emplace(vertex, 0);
    }
    std::size_t time = 0;
    std::vector<std::size_t> bridges;

    auto dfs = [&](auto &self, NetworKit::node vertex,
                   std::size_t parent_edge) -> void {
        discovery.at(vertex) = low.at(vertex) = ++time;
        for (const auto dart : graph.rotation(vertex)) {
            const auto edge = graph.dart(dart).edge;
            const auto neighbor = graph.dart(dart).to;
            if (edge == parent_edge) {
                continue;
            }
            if (discovery.at(neighbor) == 0) {
                self(self, neighbor, edge);
                low.at(vertex) = std::min(
                    low.at(vertex), low.at(neighbor));
                if (low.at(neighbor) > discovery.at(vertex)) {
                    bridges.push_back(edge);
                }
            } else {
                low.at(vertex) = std::min(
                    low.at(vertex), discovery.at(neighbor));
            }
        }
    };

    for (const auto vertex : vertices) {
        if (discovery.at(vertex) == 0) {
            dfs(dfs, vertex, invalid_index);
        }
    }
    return bridges;
}

std::vector<NetworKit::node> collectWedgeTargets(
        const ComponentWork &work, const PlaneGraph &global_graph,
        std::size_t incoming, std::size_t outgoing) {
    if (incoming >= work.referenceGlobalDart.size()
        || outgoing >= work.referenceGlobalDart.size()
        || !work.referenceGlobalDart[incoming].has_value()
        || !work.referenceGlobalDart[outgoing].has_value()) {
        return {};
    }
    const std::size_t global_incoming =
        *work.referenceGlobalDart[incoming];
    const std::size_t global_outgoing =
        *work.referenceGlobalDart[outgoing];
    std::vector<NetworKit::node> result;
    std::size_t current = global_graph.dart(
        global_graph.dart(global_incoming).twin).rotationPrevious;
    std::size_t steps = 0;
    while (current != global_outgoing) {
        result.push_back(global_graph.dart(current).to);
        current = global_graph.dart(current).rotationPrevious;
        if (++steps > global_graph.numberOfDarts()) {
            throw std::logic_error("A local face wedge did not close");
        }
    }
    return result;
}

ComponentWork makeComponentWork(
        const PlaneGraph &global_graph,
        const std::vector<NetworKit::node> &vertices,
        NetworKit::count level,
        const std::vector<std::size_t> &edge_indices,
        const std::unordered_map<
            NetworKit::node, std::size_t> &work_of_node,
        const LevelMap &levels,
        std::size_t global_outer_face,
        const PlaneGraph::Faces &global_faces,
        std::vector<Koala::BakerFakeEdge> &fake_edges) {
    ComponentWork work;
    work.level = level;
    work.vertices = vertices;
    work.embedding = PlaneGraph(global_graph.nodeBound());

    std::unordered_map<std::size_t, std::size_t> global_to_local;
    global_to_local.reserve(2 * edge_indices.size());
    for (const auto edge_index : edge_indices) {
        const auto &edge = global_graph.edge(edge_index);
        const auto global_first = edge.first;
        const auto u = global_graph.dart(global_first).from;
        const auto v = global_graph.dart(global_first).to;
        const std::size_t local_edge = work.embedding.addDetachedEdge(
            u, v, edge.original, edge.fakeType);
        const auto &created = work.embedding.edge(local_edge);
        const std::size_t local_first =
            work.embedding.dart(created.first).from == u
            ? created.first : created.second;
        const std::size_t local_second =
            work.embedding.dart(created.first).from == v
            ? created.first : created.second;
        global_to_local.emplace(global_first, local_first);
        global_to_local.emplace(
            global_graph.dart(global_first).twin, local_second);
    }

    work.globalDart.resize(work.embedding.numberOfDarts());
    work.referenceGlobalDart.resize(work.embedding.numberOfDarts());
    for (const auto &[global_dart, local_dart] : global_to_local) {
        work.globalDart[local_dart] = global_dart;
        work.referenceGlobalDart[local_dart] = global_dart;
    }
    for (const auto vertex : vertices) {
        std::vector<std::size_t> rotation;
        for (const auto global_dart : global_graph.rotation(vertex)) {
            const auto iterator = global_to_local.find(global_dart);
            if (iterator != global_to_local.end()) {
                rotation.push_back(iterator->second);
            }
        }
        work.embedding.setRotation(vertex, rotation);
    }

    const auto pre_faces = work.embedding.faces();
    std::size_t outer_face = invalid_index;
    std::vector<std::optional<std::size_t>> enclosed_work(
        pre_faces.boundaries.size());
    for (std::size_t face = 0;
         face < pre_faces.boundaries.size(); ++face) {
        std::unordered_set<std::size_t> higher_components;
        bool has_lower_neighbor = false;
        bool contains_global_outer_face = false;
        const auto &boundary = pre_faces.boundaries[face];
        for (std::size_t i = 0; i < boundary.size(); ++i) {
            const std::size_t incoming =
                boundary[(i + boundary.size() - 1) % boundary.size()];
            const std::size_t outgoing = boundary[i];
            if (work.globalDart[outgoing].has_value()
                && global_faces.faceOfDart[*work.globalDart[outgoing]]
                    == global_outer_face) {
                contains_global_outer_face = true;
            }
            for (const auto target : collectWedgeTargets(
                     work, global_graph, incoming, outgoing)) {
                if (levels.at(target) + 1 == level) {
                    has_lower_neighbor = true;
                } else if (levels.at(target) == level + 1) {
                    higher_components.insert(work_of_node.at(target));
                }
            }
        }
        if ((level == 1 && contains_global_outer_face)
            || (level > 1 && has_lower_neighbor)) {
            if (outer_face != invalid_index && outer_face != face) {
                throw std::logic_error(
                    "A level component has more than one exterior face");
            }
            outer_face = face;
        }
        if (higher_components.size() > 1) {
            throw std::logic_error(
                "Triangulation did not connect a nested level component");
        }
        if (!higher_components.empty()) {
            enclosed_work[face] = *higher_components.begin();
        }
    }
    if (outer_face == invalid_index) {
        throw std::logic_error("A level component has no exterior face");
    }

    std::size_t outer_witness = pre_faces.boundaries[outer_face].front();
    const std::vector<std::size_t> pre_face_of_dart =
        pre_faces.faceOfDart;
    for (const auto bridge : findBridges(work.embedding, vertices)) {
        const auto bridge_edge = work.embedding.edge(bridge);
        const auto bridge_dart = bridge_edge.first;
        const auto next = work.embedding.faceNext(bridge_dart);
        const auto rotation_next =
            work.embedding.dart(bridge_dart).rotationNext;
        const auto previous = work.embedding.dart(rotation_next).twin;
        const auto [copy_forward, copy_backward] =
            work.embedding.addEdgeInFace(
                previous, bridge_dart, bridge_dart, next, false,
                Koala::BakerFakeEdgeType::BRIDGE_COPY);
        const auto u = work.embedding.dart(copy_forward).from;
        const auto v = work.embedding.dart(copy_forward).to;
        const auto [edge_u, edge_v] = orderedEdge(u, v);
        fake_edges.push_back(
            {edge_u, edge_v, Koala::BakerFakeEdgeType::BRIDGE_COPY});

        work.globalDart.resize(work.embedding.numberOfDarts());
        work.referenceGlobalDart.resize(work.embedding.numberOfDarts());
        const std::size_t reference_forward =
            work.embedding.dart(bridge_edge.first).from == u
            ? bridge_edge.first : bridge_edge.second;
        if (!work.referenceGlobalDart[reference_forward].has_value()) {
            throw std::logic_error("A bridge has no global reference dart");
        }
        work.referenceGlobalDart[copy_forward] =
            work.referenceGlobalDart[reference_forward];
        work.referenceGlobalDart[copy_backward] =
            work.referenceGlobalDart[
                work.embedding.dart(reference_forward).twin];
        outer_witness = copy_forward;
    }

    const auto final_faces = work.embedding.faces();
    const std::size_t final_outer =
        final_faces.faceOfDart[outer_witness];
    std::vector<std::size_t> tree_vertex_of_face(
        final_faces.boundaries.size(), invalid_index);
    std::vector<std::optional<std::size_t>> final_enclosed_work(
        final_faces.boundaries.size());

    for (std::size_t face = 0;
         face < final_faces.boundaries.size(); ++face) {
        if (face == final_outer) {
            continue;
        }
        WorkTreeVertex tree_vertex;
        tree_vertex.face = true;
        tree_vertex.localFace = face;
        tree_vertex_of_face[face] = work.treeVertices.size();
        work.treeVertices.push_back(std::move(tree_vertex));

        std::unordered_set<std::size_t> preface_ids;
        for (const auto dart : final_faces.boundaries[face]) {
            if (dart < pre_face_of_dart.size()) {
                preface_ids.insert(pre_face_of_dart[dart]);
            }
        }
        for (const auto preface : preface_ids) {
            if (!enclosed_work[preface].has_value()) {
                continue;
            }
            if (final_enclosed_work[face].has_value()
                && final_enclosed_work[face] != enclosed_work[preface]) {
                throw std::logic_error(
                    "A face contains two nested level components");
            }
            final_enclosed_work[face] = enclosed_work[preface];
        }
    }

    std::unordered_set<
        std::pair<std::size_t, std::size_t>, IndexPairHash> tree_edges;
    auto add_tree_edge = [&](std::size_t first, std::size_t second) {
        if (first == second) {
            return;
        }
        const auto ordered = std::minmax(first, second);
        if (!tree_edges.insert(ordered).second) {
            return;
        }
        work.treeVertices[first].neighbors.push_back(second);
        work.treeVertices[second].neighbors.push_back(first);
    };

    const auto &outer_boundary = final_faces.boundaries[final_outer];
    for (std::size_t offset = 0;
         offset < outer_boundary.size(); ++offset) {
        const auto outer_dart =
            outer_boundary[outer_boundary.size() - offset - 1];
        const auto dart = work.embedding.dart(outer_dart).twin;
        const auto inner_face = final_faces.faceOfDart[dart];
        if (inner_face == final_outer
            || tree_vertex_of_face[inner_face] == invalid_index) {
            throw std::logic_error(
                "An exterior level edge has no interior face");
        }
        WorkTreeVertex leaf;
        leaf.leafDart = dart;
        const std::size_t leaf_index = work.treeVertices.size();
        work.treeVertices.push_back(std::move(leaf));
        work.outerLeaves.push_back(leaf_index);
        work.outerLeafDarts.push_back(dart);
        add_tree_edge(leaf_index, tree_vertex_of_face[inner_face]);
    }

    for (const auto &edge : work.embedding.edges()) {
        const auto first_face = final_faces.faceOfDart[edge.first];
        const auto second_face = final_faces.faceOfDart[edge.second];
        if (first_face == final_outer || second_face == final_outer
            || first_face == second_face) {
            continue;
        }
        add_tree_edge(
            tree_vertex_of_face[first_face],
            tree_vertex_of_face[second_face]);
    }

    DisjointSet tree_components(work.treeVertices.size());
    for (const auto &[first, second] : tree_edges) {
        tree_components.unite(first, second);
    }
    std::unordered_map<
        NetworKit::node, std::vector<std::size_t>> faces_at_vertex;
    faces_at_vertex.reserve(vertices.size());
    for (std::size_t face = 0;
         face < final_faces.boundaries.size(); ++face) {
        if (face == final_outer) {
            continue;
        }
        std::unordered_set<NetworKit::node> boundary_vertices;
        for (const auto dart : final_faces.boundaries[face]) {
            boundary_vertices.insert(work.embedding.dart(dart).from);
        }
        for (const auto vertex : boundary_vertices) {
            faces_at_vertex[vertex].push_back(
                tree_vertex_of_face[face]);
        }
    }
    for (const auto vertex : vertices) {
        auto candidates = std::move(faces_at_vertex[vertex]);
        if (candidates.empty()) {
            continue;
        }
        const auto first = candidates.front();
        for (std::size_t second = 1; second < candidates.size(); ++second) {
            if (tree_components.unite(first, candidates[second])) {
                add_tree_edge(first, candidates[second]);
            }
        }
    }
    if (!work.treeVertices.empty()) {
        const std::size_t representative = tree_components.find(0);
        for (std::size_t vertex = 1;
             vertex < work.treeVertices.size(); ++vertex) {
            if (tree_components.find(vertex) != representative) {
                throw std::logic_error(
                    "The ordered face construction produced a forest");
            }
        }
        if (tree_edges.size() + 1 != work.treeVertices.size()) {
            std::string details;
            for (const auto vertex : vertices) {
                details += " v" + std::to_string(vertex);
            }
            for (const auto &edge : work.embedding.edges()) {
                details += " e" + std::to_string(
                    work.embedding.dart(edge.first).from)
                    + "-" + std::to_string(
                        work.embedding.dart(edge.first).to)
                    + "[" + std::to_string(
                        final_faces.faceOfDart[edge.first])
                    + "," + std::to_string(
                        final_faces.faceOfDart[edge.second]) + "]";
            }
            details += " outer=" + std::to_string(final_outer);
            throw std::logic_error(
                "The ordered face construction is not a tree at level "
                + std::to_string(level) + ": "
                + std::to_string(work.treeVertices.size()) + " vertices, "
                + std::to_string(tree_edges.size()) + " edges; component "
                + std::to_string(vertices.size()) + " vertices, "
                + std::to_string(work.embedding.edges().size())
                + " edges, " + std::to_string(final_faces.boundaries.size())
                + " faces;" + details);
        }
    }

    for (std::size_t face = 0;
         face < final_enclosed_work.size(); ++face) {
        if (face == final_outer
            || !final_enclosed_work[face].has_value()) {
            continue;
        }
        work.treeVertices[tree_vertex_of_face[face]].enclosedWork =
            final_enclosed_work[face];
    }
    return work;
}

struct RootedWork {
    std::size_t root = invalid_index;
    std::vector<std::optional<std::size_t>> parent;
    std::vector<std::vector<std::size_t>> children;
    std::vector<NetworKit::node> label_left;
    std::vector<NetworKit::node> label_right;
    std::vector<std::size_t> ordered_leaves;
};

RootedWork rootAndOrder(
        const ComponentWork &work, std::size_t start_position) {
    if (work.outerLeaves.empty()) {
        throw std::logic_error("Cannot root an edgeless face tree");
    }
    const std::size_t leaf_count = work.outerLeaves.size();
    start_position %= leaf_count;
    const std::size_t first_leaf = work.outerLeaves[start_position];
    if (work.treeVertices[first_leaf].neighbors.size() != 1) {
        throw std::logic_error("An exterior edge vertex is not a leaf");
    }

    RootedWork result;
    result.root = work.treeVertices[first_leaf].neighbors.front();
    result.parent.resize(work.treeVertices.size());
    result.children.resize(work.treeVertices.size());
    result.label_left.assign(work.treeVertices.size(), NetworKit::none);
    result.label_right.assign(work.treeVertices.size(), NetworKit::none);
    std::vector<std::size_t> order;
    order.push_back(result.root);
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto vertex = order[i];
        for (const auto neighbor : work.treeVertices[vertex].neighbors) {
            if (result.parent[vertex].has_value()
                && neighbor == *result.parent[vertex]) {
                continue;
            }
            if (neighbor == result.root
                || result.parent[neighbor].has_value()) {
                throw std::logic_error("A face tree contains a cycle");
            }
            result.parent[neighbor] = vertex;
            order.push_back(neighbor);
        }
    }

    std::vector<std::size_t> rotated_position(
        work.treeVertices.size(), invalid_index);
    result.ordered_leaves.reserve(leaf_count);
    for (std::size_t position = 0; position < leaf_count; ++position) {
        const std::size_t original =
            (start_position + position) % leaf_count;
        const auto leaf = work.outerLeaves[original];
        rotated_position[leaf] = position;
        result.ordered_leaves.push_back(leaf);
    }

    std::vector<std::size_t> minimum(
        work.treeVertices.size(), invalid_index);
    std::vector<std::size_t> maximum(work.treeVertices.size(), 0);
    std::vector<std::size_t> count(work.treeVertices.size(), 0);
    for (auto iterator = order.rbegin(); iterator != order.rend(); ++iterator) {
        const auto vertex = *iterator;
        if (!work.treeVertices[vertex].face) {
            minimum[vertex] = maximum[vertex] = rotated_position[vertex];
            count[vertex] = 1;
        }
        for (const auto neighbor : work.treeVertices[vertex].neighbors) {
            if (!result.parent[neighbor].has_value()
                || *result.parent[neighbor] != vertex) {
                continue;
            }
            if (count[neighbor] == 0) {
                throw std::logic_error(
                    "A face-tree subtree contains no exterior edge");
            }
            result.children[vertex].push_back(neighbor);
            minimum[vertex] = std::min(minimum[vertex], minimum[neighbor]);
            maximum[vertex] = std::max(maximum[vertex], maximum[neighbor]);
            count[vertex] += count[neighbor];
        }
        if (!result.children[vertex].empty()) {
            std::unordered_map<std::size_t, std::size_t> child_at_position;
            child_at_position.reserve(result.children[vertex].size());
            for (const auto child : result.children[vertex]) {
                if (!child_at_position.emplace(
                        minimum[child], child).second) {
                    throw std::logic_error(
                        "Two face-tree children start at one leaf");
                }
            }
            std::vector<std::size_t> ordered_children;
            ordered_children.reserve(result.children[vertex].size());
            std::size_t position = minimum[vertex];
            while (ordered_children.size()
                   < result.children[vertex].size()) {
                const auto child = child_at_position.find(position);
                if (child == child_at_position.end()) {
                    throw std::logic_error(
                        "Face-tree child intervals contain a gap");
                }
                ordered_children.push_back(child->second);
                position = maximum[child->second] + 1;
            }
            result.children[vertex] = std::move(ordered_children);
        }
        if (count[vertex] != 0
            && maximum[vertex] - minimum[vertex] + 1 != count[vertex]) {
            throw std::logic_error(
                "A face-tree subtree is not an exterior-edge interval");
        }
    }

    for (auto iterator = order.rbegin(); iterator != order.rend(); ++iterator) {
        const auto vertex = *iterator;
        if (!work.treeVertices[vertex].face) {
            const auto dart = work.treeVertices[vertex].leafDart;
            result.label_left[vertex] = work.embedding.dart(dart).from;
            result.label_right[vertex] = work.embedding.dart(dart).to;
            continue;
        }
        if (result.children[vertex].empty()) {
            throw std::logic_error("A rooted face vertex has no children");
        }
        result.label_left[vertex] =
            result.label_left[result.children[vertex].front()];
        result.label_right[vertex] =
            result.label_right[result.children[vertex].back()];
        for (std::size_t i = 1;
             i < result.children[vertex].size(); ++i) {
            const auto previous = result.children[vertex][i - 1];
            const auto current = result.children[vertex][i];
            if (result.label_right[previous]
                != result.label_left[current]) {
                throw std::logic_error(
                    "Ordered face-tree child labels do not form a walk");
            }
        }
    }
    if (result.label_left[result.root] != result.label_right[result.root]) {
        throw std::logic_error("An ordered face-tree root is not labelled (z,z)");
    }
    return result;
}

std::size_t selectStartPosition(
        const ComponentWork &work, const PlaneGraph &global_graph,
        NetworKit::node parent_left, NetworKit::node parent_right) {
    NetworKit::node root_vertex = NetworKit::none;
    std::optional<std::size_t> cross_dart;
    if (parent_left != parent_right) {
        for (const auto candidate : work.vertices) {
            if (!global_graph.hasEdge(candidate, parent_left)
                || !global_graph.hasEdge(candidate, parent_right)) {
                continue;
            }
            for (const auto candidate_dart :
                 global_graph.rotation(candidate)) {
                if (global_graph.dart(candidate_dart).to != parent_left) {
                    continue;
                }
                for (const auto face_dart : {
                         candidate_dart,
                         global_graph.dart(candidate_dart).twin}) {
                    const auto second_dart =
                        global_graph.faceNext(face_dart);
                    if (global_graph.dart(second_dart).to
                        != parent_right) {
                        continue;
                    }
                    const auto third_dart =
                        global_graph.faceNext(second_dart);
                    if (global_graph.dart(third_dart).to
                            == global_graph.dart(face_dart).from
                        && global_graph.faceNext(third_dart)
                            == face_dart) {
                        root_vertex = candidate;
                        cross_dart = candidate_dart;
                        break;
                    }
                }
                if (root_vertex != NetworKit::none) {
                    break;
                }
            }
            if (root_vertex != NetworKit::none) {
                break;
            }
        }
    } else {
        for (const auto vertex : work.vertices) {
            if (global_graph.hasEdge(vertex, parent_left)) {
                root_vertex = std::min(root_vertex, vertex);
            }
        }
    }
    if (root_vertex == NetworKit::none) {
        throw std::logic_error(
            "No level-tree root is adjacent to its enclosing label");
    }

    std::unordered_map<std::size_t, std::size_t> position_of_reference;
    for (std::size_t position = 0;
         position < work.outerLeafDarts.size(); ++position) {
        const auto local_dart = work.outerLeafDarts[position];
        if (work.embedding.dart(local_dart).from != root_vertex
            || local_dart >= work.referenceGlobalDart.size()
            || !work.referenceGlobalDart[local_dart].has_value()) {
            continue;
        }
        position_of_reference.emplace(
            *work.referenceGlobalDart[local_dart], position);
    }
    if (position_of_reference.empty()) {
        throw std::logic_error(
            "A level-tree root has no oriented exterior edge");
    }

    if (!cross_dart.has_value()) {
        for (const auto dart : global_graph.rotation(root_vertex)) {
            if (global_graph.dart(dart).to == parent_left) {
                cross_dart = dart;
                break;
            }
        }
    }
    if (!cross_dart.has_value()) {
        throw std::logic_error(
            "A level-tree root has no edge to its enclosing label");
    }
    std::size_t current =
        global_graph.dart(*cross_dart).rotationNext;
    while (current != *cross_dart) {
        const auto iterator = position_of_reference.find(current);
        if (iterator != position_of_reference.end()) {
            return iterator->second;
        }
        current = global_graph.dart(current).rotationNext;
    }
    throw std::logic_error(
        "No exterior edge follows the enclosing cross edge");
}

std::vector<NetworKit::node> boundaryAt(
        const Koala::BakerFaceTree &tree,
        const Koala::BakerFaceTreeVertex &face,
        std::size_t boundary_index) {
    if (face.children.empty()
        || boundary_index > face.children.size()) {
        throw std::logic_error("Invalid LB/RB boundary index");
    }
    if (boundary_index == face.children.size()) {
        return tree.vertices[face.children.back()].boundary.right;
    }
    return tree.vertices[face.children[boundary_index]].boundary.left;
}

void verifyBoundaryLevels(
        const Koala::BakerFaceTree &tree,
        const Koala::BakerFaceTreeVertex &vertex,
        const std::vector<NetworKit::count> &levels) {
    if (vertex.boundary.left.size() != tree.level
        || vertex.boundary.right.size() != tree.level) {
        throw std::logic_error("A Baker boundary does not have length i");
    }
    for (std::size_t i = 0; i < tree.level; ++i) {
        const auto expected = tree.level - i;
        if (levels[vertex.boundary.left[i]] != expected
            || levels[vertex.boundary.right[i]] != expected) {
            throw std::logic_error(
                "A Baker boundary does not contain one vertex per level");
        }
    }
}

class ForestBuilder {
 public:
    ForestBuilder(
        const NetworKit::Graph &input, PlaneGraph global_graph,
        LevelMap levels,
        std::size_t outer_face, PlaneGraph::Faces global_faces,
        std::vector<NetworKit::node> component,
        Koala::BakerForest &forest)
        : input_(input),
          globalGraph_(std::move(global_graph)),
          levels_(std::move(levels)),
          outerFace_(outer_face),
          globalFaces_(std::move(global_faces)),
          component_(std::move(component)),
          forest_(forest) { }

    void build() {
        identifyLevelComponents();
        buildComponentWorks();
        identifyParents();
        std::vector<std::size_t> top_components;
        for (std::size_t work = 0; work < works_.size(); ++work) {
            if (works_[work].level == 1) {
                top_components.push_back(work);
            }
        }
        if (top_components.size() != 1) {
            throw std::logic_error(
                "A connected embedding has more than one level-1 component");
        }
        const std::size_t root_tree =
            materialize(top_components.front(), std::nullopt, std::nullopt);
        forest_.roots.push_back(root_tree);
    }

 private:
    void identifyLevelComponents() {
        DisjointSet components(globalGraph_.nodeBound());
        for (const auto &edge : globalGraph_.edges()) {
            const auto u = globalGraph_.dart(edge.first).from;
            const auto v = globalGraph_.dart(edge.first).to;
            if (levels_.at(u) == levels_.at(v)) {
                components.unite(u, v);
            }
        }

        std::unordered_map<std::size_t, std::size_t> ids;
        workOfNode_.clear();
        workOfNode_.reserve(component_.size());
        for (const auto vertex : component_) {
            const auto key = components.find(vertex);
            auto [iterator, inserted] = ids.emplace(key, ids.size());
            static_cast<void>(inserted);
            workOfNode_.emplace(vertex, iterator->second);
        }
        workVertices_.resize(ids.size());
        workLevels_.resize(ids.size());
        workEdges_.resize(ids.size());
        for (const auto vertex : component_) {
            const auto work = workOfNode_.at(vertex);
            workVertices_[work].push_back(vertex);
            workLevels_[work] = levels_.at(vertex);
        }
        for (std::size_t edge = 0;
             edge < globalGraph_.edges().size(); ++edge) {
            const auto &embedded_edge = globalGraph_.edge(edge);
            const auto u = globalGraph_.dart(embedded_edge.first).from;
            const auto v = globalGraph_.dart(embedded_edge.first).to;
            if (levels_.at(u) == levels_.at(v)
                && workOfNode_.at(u) == workOfNode_.at(v)) {
                workEdges_[workOfNode_.at(u)].push_back(edge);
            }
        }
    }

    void buildComponentWorks() {
        works_.resize(workVertices_.size());
        for (std::size_t work = 0; work < works_.size(); ++work) {
            if (workVertices_[work].size() == 1
                && !hasSameLevelEdge(workVertices_[work].front())) {
                works_[work].level = workLevels_[work];
                works_[work].vertices = workVertices_[work];
                continue;
            }
            works_[work] = makeComponentWork(
                globalGraph_, workVertices_[work], workLevels_[work],
                workEdges_[work], workOfNode_, levels_, outerFace_, globalFaces_,
                forest_.fakeEdges);
        }
    }

    bool hasSameLevelEdge(NetworKit::node vertex) const {
        for (const auto dart : globalGraph_.rotation(vertex)) {
            if (levels_.at(globalGraph_.dart(dart).to)
                == levels_.at(vertex)) {
                return true;
            }
        }
        return false;
    }

    void identifyParents() {
        for (std::size_t parent = 0; parent < works_.size(); ++parent) {
            for (std::size_t vertex = 0;
                 vertex < works_[parent].treeVertices.size(); ++vertex) {
                const auto child =
                    works_[parent].treeVertices[vertex].enclosedWork;
                if (!child.has_value()) {
                    continue;
                }
                if (works_[*child].parentWork.has_value()) {
                    throw std::logic_error(
                        "A level component has two enclosing faces");
                }
                works_[*child].parentWork = parent;
                works_[*child].parentFaceVertex = vertex;
            }
        }

        for (std::size_t child = 0; child < works_.size(); ++child) {
            if (works_[child].level == 1) {
                continue;
            }
            if (works_[child].parentWork.has_value()) {
                continue;
            }
            const auto vertex = works_[child].vertices.front();
            std::unordered_set<std::size_t> possible_parents;
            for (const auto dart : globalGraph_.rotation(vertex)) {
                const auto neighbor = globalGraph_.dart(dart).to;
                if (levels_.at(neighbor) + 1
                    == levels_.at(vertex)) {
                    possible_parents.insert(workOfNode_.at(neighbor));
                }
            }
            if (possible_parents.size() != 1) {
                throw std::logic_error(
                    "A nested singleton has no unique enclosing component");
            }
            const auto parent = *possible_parents.begin();
            std::optional<std::size_t> parent_face;
            const auto parent_faces = works_[parent].embedding.faces();
            for (std::size_t index = 0;
                 index < works_[parent].treeVertices.size(); ++index) {
                if (!works_[parent].treeVertices[index].face) {
                    continue;
                }
                const auto local_face =
                    works_[parent].treeVertices[index].localFace;
                bool incident = false;
                for (const auto boundary_dart :
                     parent_faces.boundaries[local_face]) {
                    const auto boundary_vertex =
                        works_[parent].embedding.dart(boundary_dart).from;
                    if (globalGraph_.hasEdge(boundary_vertex, vertex)) {
                        incident = true;
                    }
                }
                if (incident) {
                    parent_face = index;
                    break;
                }
            }
            if (!parent_face.has_value()) {
                throw std::logic_error(
                    "A nested singleton has no enclosing face vertex");
            }
            works_[child].parentWork = parent;
            works_[child].parentFaceVertex = parent_face;
            works_[parent].treeVertices[*parent_face].enclosedWork = child;
        }
    }

    std::size_t materialize(
            std::size_t work_index,
            std::optional<std::size_t> enclosing_tree,
            std::optional<std::size_t> enclosing_vertex) {
        auto &work = works_[work_index];
        if (work.publicTree.has_value()) {
            throw std::logic_error("A Baker level tree was materialized twice");
        }

        const bool singleton = work.treeVertices.empty();
        RootedWork rooted;
        if (!singleton) {
            std::size_t start = 0;
            if (enclosing_tree.has_value()) {
                const auto &parent_vertex =
                    forest_.trees[*enclosing_tree].vertices[*enclosing_vertex];
                start = selectStartPosition(
                    work, globalGraph_, parent_vertex.labelLeft,
                    parent_vertex.labelRight);
            }
            rooted = rootAndOrder(work, start);
        }

        Koala::BakerFaceTree tree;
        tree.level = work.level;
        tree.component = work.vertices;
        tree.enclosingTree = enclosing_tree;
        tree.enclosingVertex = enclosing_vertex;
        if (singleton) {
            Koala::BakerFaceTreeVertex root;
            root.type = Koala::BakerFaceTreeVertexType::SINGLETON;
            root.labelLeft = work.vertices.front();
            root.labelRight = work.vertices.front();
            tree.vertices.push_back(std::move(root));
            tree.root = 0;
        } else {
            tree.vertices.resize(work.treeVertices.size());
            tree.root = rooted.root;
            tree.leaves = rooted.ordered_leaves;
            for (std::size_t vertex = 0;
                 vertex < work.treeVertices.size(); ++vertex) {
                auto &output = tree.vertices[vertex];
                output.type = work.treeVertices[vertex].face
                    ? Koala::BakerFaceTreeVertexType::FACE
                    : Koala::BakerFaceTreeVertexType::EXTERIOR_EDGE;
                output.labelLeft = rooted.label_left[vertex];
                output.labelRight = rooted.label_right[vertex];
                output.parent = rooted.parent[vertex];
                output.children = rooted.children[vertex];
                if (!work.treeVertices[vertex].face) {
                    const auto local_edge = work.embedding.dart(
                        work.treeVertices[vertex].leafDart).edge;
                    output.originalLabelEdge =
                        work.embedding.edge(local_edge).original;
                } else {
                    output.originalLabelEdge =
                        output.labelLeft != output.labelRight
                        && input_.hasEdge(
                            output.labelLeft, output.labelRight);
                }
            }
        }

        const std::size_t public_tree = forest_.trees.size();
        forest_.trees.push_back(std::move(tree));
        work.publicTree = public_tree;
        computeBoundaries(public_tree, work);

        if (!singleton) {
            for (std::size_t vertex = 0;
                 vertex < work.treeVertices.size(); ++vertex) {
                const auto child =
                    work.treeVertices[vertex].enclosedWork;
                if (!child.has_value()) {
                    continue;
                }
                const std::size_t child_tree =
                    materialize(*child, public_tree, vertex);
                forest_.trees[public_tree].vertices[vertex].enclosedTree =
                    child_tree;
            }
        }
        return public_tree;
    }

    void computeBoundaries(
            std::size_t tree_index, const ComponentWork &work) {
        auto &tree = forest_.trees[tree_index];
        if (tree.level == 1) {
            computeLevelOneBoundaries(tree);
            return;
        }
        if (!tree.enclosingTree.has_value()
            || !tree.enclosingVertex.has_value()) {
            throw std::logic_error("A deeper Baker tree has no enclosing face");
        }
        const auto &parent_tree = forest_.trees[*tree.enclosingTree];
        const auto &parent_face =
            parent_tree.vertices[*tree.enclosingVertex];
        const std::size_t number_of_cuts = parent_face.children.size();
        if (number_of_cuts == 0) {
            throw std::logic_error("An enclosing face has no ordered children");
        }

        if (tree.vertices[tree.root].type
            == Koala::BakerFaceTreeVertexType::SINGLETON) {
            auto &root = tree.vertices[tree.root];
            root.leftBoundaryIndex = 0;
            root.rightBoundaryIndex = number_of_cuts;
        } else {
            computeLeafBoundaryIndices(
                tree, parent_tree, parent_face, work);
            std::vector<std::size_t> order;
            order.push_back(tree.root);
            for (std::size_t i = 0; i < order.size(); ++i) {
                for (const auto child : tree.vertices[order[i]].children) {
                    order.push_back(child);
                }
            }
            for (auto iterator = order.rbegin();
                 iterator != order.rend(); ++iterator) {
                auto &vertex = tree.vertices[*iterator];
                if (vertex.type
                    != Koala::BakerFaceTreeVertexType::FACE) {
                    continue;
                }
                vertex.leftBoundaryIndex =
                    tree.vertices[vertex.children.front()].leftBoundaryIndex;
                vertex.rightBoundaryIndex =
                    tree.vertices[vertex.children.back()].rightBoundaryIndex;
            }
        }
        computeMiddleBoundaryIndices(
            tree, parent_tree, parent_face, work);

        for (auto &vertex : tree.vertices) {
            vertex.boundary.left.push_back(vertex.labelLeft);
            vertex.boundary.right.push_back(vertex.labelRight);
            auto left_suffix = boundaryAt(
                parent_tree, parent_face, vertex.leftBoundaryIndex);
            auto right_suffix = boundaryAt(
                parent_tree, parent_face, vertex.rightBoundaryIndex);
            vertex.boundary.left.insert(
                vertex.boundary.left.end(),
                left_suffix.begin(), left_suffix.end());
            vertex.boundary.right.insert(
                vertex.boundary.right.end(),
                right_suffix.begin(), right_suffix.end());
            verifyBoundaryLevels(tree, vertex, forest_.levels);
        }
        verifySiblingBoundaries(tree);
    }

    NetworKit::node lowerFaceVertex(
            const ComponentWork &work, std::size_t vertex_index,
            NetworKit::count level) const {
        const auto leaf_dart =
            work.treeVertices[vertex_index].leafDart;
        const auto outer_dart =
            work.embedding.dart(leaf_dart).twin;
        if (outer_dart >= work.referenceGlobalDart.size()
            || !work.referenceGlobalDart[outer_dart].has_value()) {
            throw std::logic_error(
                "A Baker exterior edge has no triangulation side");
        }
        const auto global_outer =
            *work.referenceGlobalDart[outer_dart];
        const auto face = globalFaces_.faceOfDart[global_outer];
        NetworKit::node result = NetworKit::none;
        for (const auto dart : globalFaces_.boundaries[face]) {
            const auto candidate = globalGraph_.dart(dart).from;
            if (levels_.at(candidate) + 1 != level) {
                continue;
            }
            if (result != NetworKit::none && result != candidate) {
                throw std::logic_error(
                    "A Baker exterior edge has two lower face vertices");
            }
            result = candidate;
        }
        if (result == NetworKit::none) {
            throw std::logic_error(
                "A Baker exterior edge has no lower face vertex");
        }
        return result;
    }

    void computeMiddleBoundaryIndices(
            Koala::BakerFaceTree &tree,
            const Koala::BakerFaceTree &parent_tree,
            const Koala::BakerFaceTreeVertex &parent_face,
            const ComponentWork &work) {
        std::vector<NetworKit::node> cuts;
        cuts.reserve(parent_face.children.size() + 1);
        for (const auto child : parent_face.children) {
            cuts.push_back(parent_tree.vertices[child].labelLeft);
        }
        cuts.push_back(
            parent_tree.vertices[parent_face.children.back()].labelRight);

        for (std::size_t vertex_index = 0;
             vertex_index < tree.vertices.size(); ++vertex_index) {
            auto &vertex = tree.vertices[vertex_index];
            if (vertex.type == Koala::BakerFaceTreeVertexType::FACE) {
                continue;
            }
            const auto left = vertex.leftBoundaryIndex;
            const auto right = vertex.rightBoundaryIndex;
            if (left > right || right >= cuts.size()) {
                throw std::logic_error("Invalid Baker LB/RB interval");
            }
            if (left == right) {
                vertex.middleBoundaryIndex = left;
                continue;
            }
            if (vertex.type
                == Koala::BakerFaceTreeVertexType::SINGLETON) {
                std::size_t selected = right;
                for (std::size_t index = left;
                     index <= right; ++index) {
                    if (globalGraph_.hasEdge(
                            vertex.labelRight, cuts[index])) {
                        selected = index;
                        break;
                    }
                }
                vertex.middleBoundaryIndex = selected;
                continue;
            }
            const auto middle_point =
                lowerFaceVertex(work, vertex_index, tree.level);

            // The incident triangle identifies Baker's occurrence of z_p.
            std::size_t selected = right;
            bool found = false;
            for (std::size_t index = left;
                 index <= right; ++index) {
                if (cuts[index] == middle_point) {
                    selected = index;
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::logic_error(
                    "A Baker middle point lies outside its LB/RB interval");
            }
            vertex.middleBoundaryIndex = selected;
        }
    }

    void computeLevelOneBoundaries(Koala::BakerFaceTree &tree) {
        std::vector<std::size_t> order;
        order.push_back(tree.root);
        for (std::size_t i = 0; i < order.size(); ++i) {
            for (const auto child : tree.vertices[order[i]].children) {
                order.push_back(child);
            }
        }
        for (auto iterator = order.rbegin();
             iterator != order.rend(); ++iterator) {
            auto &vertex = tree.vertices[*iterator];
            if (vertex.type == Koala::BakerFaceTreeVertexType::FACE) {
                vertex.boundary.left =
                    tree.vertices[vertex.children.front()].boundary.left;
                vertex.boundary.right =
                    tree.vertices[vertex.children.back()].boundary.right;
            } else {
                vertex.boundary.left = {vertex.labelLeft};
                vertex.boundary.right = {vertex.labelRight};
            }
            verifyBoundaryLevels(tree, vertex, forest_.levels);
        }
        verifySiblingBoundaries(tree);
    }

    void computeLeafBoundaryIndices(
            Koala::BakerFaceTree &tree,
            const Koala::BakerFaceTree &parent_tree,
            const Koala::BakerFaceTreeVertex &parent_face,
            const ComponentWork &work) {
        std::vector<NetworKit::node> cuts;
        cuts.reserve(parent_face.children.size() + 1);
        for (const auto child : parent_face.children) {
            cuts.push_back(parent_tree.vertices[child].labelLeft);
        }
        cuts.push_back(
            parent_tree.vertices[parent_face.children.back()].labelRight);
        std::unordered_map<
            NetworKit::node, std::vector<std::size_t>> cut_positions;
        for (std::size_t cut = 0; cut < cuts.size(); ++cut) {
            cut_positions[cuts[cut]].push_back(cut);
        }

        if (tree.leaves.empty()) {
            throw std::logic_error("A non-singleton Baker tree has no leaves");
        }
        std::vector<std::vector<NetworKit::node>> transition_points(
            tree.leaves.size());
        for (std::size_t leaf = 1; leaf < tree.leaves.size(); ++leaf) {
            const auto previous = tree.leaves[leaf - 1];
            const auto current = tree.leaves[leaf];
            const auto incoming = work.embedding.dart(
                work.treeVertices[current].leafDart).twin;
            const auto outgoing = work.embedding.dart(
                work.treeVertices[previous].leafDart).twin;
            transition_points[leaf] = collectWedgeTargets(
                work, globalGraph_, incoming, outgoing);
        }

        std::unordered_map<NetworKit::node, std::size_t> cut_cursors;
        std::unordered_map<NetworKit::node, std::size_t> original_cursors;
        tree.vertices[tree.leaves.front()].leftBoundaryIndex = 0;
        for (std::size_t leaf = 1; leaf < tree.leaves.size(); ++leaf) {
            const auto previous = tree.leaves[leaf - 1];
            const auto current = tree.leaves[leaf];
            const auto dividing_vertex =
                tree.vertices[current].labelLeft;
            if (tree.vertices[previous].labelRight != dividing_vertex) {
                throw std::logic_error(
                    "Successive level-tree leaves do not share a node");
            }
            const auto &dividing_points = transition_points[leaf];
            const std::size_t lower_bound =
                tree.vertices[previous].leftBoundaryIndex;
            std::size_t selected = invalid_index;
            // The plane graph can contain several occurrences of one cutpoint.
            // A vertex-only lookup cannot distinguish Baker's y_p occurrences.
            // Prefer the next distinct original dividing point when one exists;
            // this advances to its actual occurrence instead of reusing y_p at
            // the left boundary through a parallel triangulation or bridge dart.
            for (const auto point : dividing_points) {
                if (point == cuts[lower_bound]
                    || !input_.hasEdge(dividing_vertex, point)) {
                    continue;
                }
                const auto positions = cut_positions.find(point);
                if (positions == cut_positions.end()) {
                    continue;
                }
                auto &cursor = original_cursors[point];
                while (cursor < positions->second.size()
                       && positions->second[cursor] <= lower_bound) {
                    ++cursor;
                }
                if (cursor < positions->second.size()) {
                    selected = std::min(
                        selected, positions->second[cursor]);
                }
            }
            for (const auto point : dividing_points) {
                if (selected != invalid_index) {
                    break;
                }
                const auto positions = cut_positions.find(point);
                if (positions == cut_positions.end()) {
                    continue;
                }
                auto &cursor = cut_cursors[point];
                while (cursor < positions->second.size()
                       && positions->second[cursor] < lower_bound) {
                    ++cursor;
                }
                if (cursor < positions->second.size()) {
                    selected = std::min(
                        selected, positions->second[cursor]);
                }
            }
            if (selected == invalid_index) {
                std::string details =
                    " at level " + std::to_string(tree.level)
                    + " between " + std::to_string(
                        tree.vertices[previous].labelLeft)
                    + "-" + std::to_string(dividing_vertex)
                    + "-" + std::to_string(
                        tree.vertices[current].labelRight)
                    + "; cuts";
                for (const auto cut : cuts) {
                    details += " " + std::to_string(cut);
                }
                details += "; wedge";
                for (const auto point : dividing_points) {
                    details += " " + std::to_string(point);
                }
                throw std::logic_error(
                    "No dividing point exists in the triangulated annulus"
                    + details);
            }
            tree.vertices[current].leftBoundaryIndex = selected;
            tree.vertices[previous].rightBoundaryIndex = selected;
        }
        tree.vertices[tree.leaves.back()].rightBoundaryIndex =
            parent_face.children.size();
    }

    void verifySiblingBoundaries(
            const Koala::BakerFaceTree &tree) const {
        for (const auto &vertex : tree.vertices) {
            for (std::size_t i = 1; i < vertex.children.size(); ++i) {
                const auto &previous =
                    tree.vertices[vertex.children[i - 1]];
                const auto &current =
                    tree.vertices[vertex.children[i]];
                if (previous.boundary.right != current.boundary.left) {
                    throw std::logic_error(
                        "Successive Baker slices have different boundaries");
                }
            }
        }
    }

    const NetworKit::Graph &input_;
    PlaneGraph globalGraph_;
    LevelMap levels_;
    std::size_t outerFace_;
    PlaneGraph::Faces globalFaces_;
    std::vector<NetworKit::node> component_;
    Koala::BakerForest &forest_;
    std::unordered_map<NetworKit::node, std::size_t> workOfNode_;
    std::vector<std::vector<NetworKit::node>> workVertices_;
    std::vector<NetworKit::count> workLevels_;
    std::vector<std::vector<std::size_t>> workEdges_;
    std::vector<ComponentWork> works_;
};

}  // namespace

namespace Koala {

BakerForest buildBakerForest(const NetworKit::Graph &graph) {
    return buildBakerForest(graph, BakerPlaneEmbedding{}, {});
}

BakerForest buildBakerForest(
        const NetworKit::Graph &graph,
        const std::vector<NetworKit::node> &outer_face_vertices) {
    return buildBakerForest(
        graph, BakerPlaneEmbedding{}, outer_face_vertices);
}

BakerForest buildBakerForest(
        const NetworKit::Graph &graph,
        const BakerPlaneEmbedding &prescribed_embedding,
        const std::vector<NetworKit::node> &outer_face_vertices) {
    if (graph.isDirected()) {
        throw std::invalid_argument(
            "BakerKOuterplanarGraphScheme requires an undirected graph");
    }

    BakerForest forest;
    forest.levels.assign(graph.upperNodeIdBound(), 0);
    forest.embedding.rotation.resize(graph.upperNodeIdBound());
    std::unordered_set<NetworKit::node> requested_outer_vertices;
    requested_outer_vertices.reserve(outer_face_vertices.size());
    for (const auto vertex : outer_face_vertices) {
        if (!graph.hasNode(vertex)) {
            throw std::invalid_argument(
                "A requested Baker outer-face vertex is missing");
        }
        requested_outer_vertices.insert(vertex);
    }
    if (graph.numberOfNodes() == 0) {
        return forest;
    }

    NetworKit::ConnectedComponents connected_components(graph);
    connected_components.run();
    for (auto component : connected_components.getComponents()) {
        std::vector<NetworKit::node> component_outer_face;
        if (!requested_outer_vertices.empty()) {
            for (const auto vertex : component) {
                if (requested_outer_vertices.contains(vertex)) {
                    component_outer_face.push_back(vertex);
                }
            }
            if (component_outer_face.empty()) {
                throw std::invalid_argument(
                    "A graph component has no requested Baker outer face");
            }
        }
        if (component.size() == 1) {
            const auto vertex = component.front();
            forest.levels[vertex] = 1;
            forest.outerplanarity_ = std::max<NetworKit::count>(
                forest.outerplanarity_, 1);
            BakerFaceTree tree;
            tree.level = 1;
            tree.component = component;
            tree.root = 0;
            BakerFaceTreeVertex root;
            root.type = BakerFaceTreeVertexType::SINGLETON;
            root.labelLeft = vertex;
            root.labelRight = vertex;
            root.boundary.left = {vertex};
            root.boundary.right = {vertex};
            tree.vertices.push_back(std::move(root));
            forest.roots.push_back(forest.trees.size());
            forest.trees.push_back(std::move(tree));
            continue;
        }

        std::vector<NetworKit::node> local_to_global;
        BoostGraph boost_graph =
            makeBoostGraph(graph, component, local_to_global);
        Embedding boost_embedding(num_vertices(boost_graph));
        const bool planar = boost::boyer_myrvold_planarity_test(
            boost::boyer_myrvold_params::graph = boost_graph,
            boost::boyer_myrvold_params::embedding =
                boost_embedding.data());
        if (!planar) {
            throw std::invalid_argument(
                "BakerKOuterplanarGraphScheme requires a planar graph");
        }
        if (!prescribed_embedding.rotation.empty()) {
            boost_embedding = restrictEmbedding(
                boost_graph, local_to_global, prescribed_embedding);
        }

        PlaneGraph plane_graph = makePlaneGraph(
            boost_graph, boost_embedding, local_to_global,
            graph.upperNodeIdBound());
        const auto initial_faces = plane_graph.faces();
        const std::size_t outer_face = component_outer_face.empty()
            ? selectOuterFace(plane_graph, initial_faces)
            : selectOuterFace(
                plane_graph, initial_faces, component_outer_face);
        auto levels = computeLevels(
            plane_graph, initial_faces, outer_face);
        for (const auto vertex : component) {
            for (const auto dart : plane_graph.rotation(vertex)) {
                forest.embedding.rotation[vertex].push_back(
                    plane_graph.dart(dart).to);
            }
        }
        for (const auto vertex : component) {
            const auto level = levels.find(vertex);
            if (level == levels.end() || level->second == 0) {
                throw std::logic_error("A planar vertex received no Baker level");
            }
            forest.levels[vertex] = level->second;
            forest.outerplanarity_ = std::max(
                forest.outerplanarity_, level->second);
        }

        const std::size_t outer_witness =
            initial_faces.boundaries[outer_face].front();
        triangulateInternalFaces(
            plane_graph, outer_face, initial_faces, levels,
            forest.fakeEdges);
        const auto triangulated_faces = plane_graph.faces();
        const std::size_t triangulated_outer =
            triangulated_faces.faceOfDart[outer_witness];

        ForestBuilder builder(
            graph, std::move(plane_graph), std::move(levels),
            triangulated_outer, triangulated_faces,
            std::move(component), forest);
        builder.build();
    }

    if (!forest.hasTwoKBoundaryBound()) {
        throw std::logic_error(
            "The materialized Baker forest violates the 2k boundary");
    }
    return forest;
}

}  // namespace Koala
