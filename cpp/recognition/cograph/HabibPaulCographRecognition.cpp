#include <algorithm>
#include <deque>
#include <list>
#include <tuple>
#include <vector>

#include "recognition/CographRecognition.hpp"

namespace Koala {

namespace {

struct PartitionPart {
    std::list<NetworKit::node> vertices;
    NetworKit::node pivot = NetworKit::none;
    NetworKit::count amount = 0;
    NetworKit::node division = NetworKit::none;
    bool active = true;
};

enum class TwinType : NetworKit::count {
    FALSE_TWINS = 0,
    TRUE_TWINS = 1,
    NOT_TWINS = 2
};

using TwinOrder = std::tuple<NetworKit::node, NetworKit::node, NodeType>;

struct PartitionRefinement {
    std::vector<PartitionPart> parts;
    std::list<NetworKit::node> order;
    std::vector<std::list<NetworKit::node>::iterator> part_position;
    std::vector<NetworKit::node> vertex_part;
    std::vector<std::list<NetworKit::node>::iterator> vertex_position;
    std::deque<NetworKit::node> unused_parts;
    NetworKit::node left = NetworKit::none;
    NetworKit::node right = NetworKit::none;
    NetworKit::node origin = NetworKit::none;
    NetworKit::count active_parts = 0;

    explicit PartitionRefinement(const NetworKit::Graph &graph)
            : vertex_part(graph.upperNodeIdBound(), NetworKit::none),
              vertex_position(graph.upperNodeIdBound()) {
        parts.emplace_back();
        part_position.resize(1);
        order.push_back(0);
        part_position[0] = order.begin();
        for (auto u : graph.nodeRange()) {
            add_vertex(0, u);
            if (origin == NetworKit::none) {
                origin = u;
            }
        }
        active_parts = graph.numberOfNodes() == 0 ? 0 : 1;
        left = right = 0;
    }

    NetworKit::node add_part_after(NetworKit::node part) {
        NetworKit::node new_part = parts.size();
        parts.emplace_back();
        part_position.push_back(order.insert(std::next(part_position[part]), new_part));
        return new_part;
    }

    void add_vertex(NetworKit::node part, NetworKit::node v) {
        parts[part].vertices.push_back(v);
        vertex_part[v] = part;
        vertex_position[v] = std::prev(parts[part].vertices.end());
    }

    void move_vertex(NetworKit::node v, NetworKit::node part) {
        auto old_part = vertex_part[v];
        parts[old_part].vertices.erase(vertex_position[v]);
        if (parts[old_part].pivot == v) {
            parts[old_part].pivot = NetworKit::none;
        }
        add_vertex(part, v);
    }

    void erase_part(NetworKit::node part) {
        if (!parts[part].active) {
            return;
        }
        order.erase(part_position[part]);
        parts[part].active = false;
        parts[part].division = NetworKit::none;
        parts[part].amount = 0;
    }

    NetworKit::node previous_part(NetworKit::node part) const {
        if (part == NetworKit::none || part_position[part] == order.begin()) {
            return NetworKit::none;
        }
        return *std::prev(part_position[part]);
    }

    NetworKit::node next_part(NetworKit::node part) const {
        if (part == NetworKit::none) {
            return NetworKit::none;
        }
        auto next = std::next(part_position[part]);
        return next == order.end() ? NetworKit::none : *next;
    }

    void check_bound(NetworKit::node &bound, bool left_side) {
        auto origin_part = vertex_part[origin];
        auto candidate = left_side ? next_part(bound) : previous_part(bound);
        if (origin_part != bound && candidate != NetworKit::none
                && parts[candidate].vertices.size() > 1) {
            bound = candidate;
        }
        candidate = left_side ? previous_part(bound) : next_part(bound);
        while (candidate != NetworKit::none && parts[bound].vertices.size() <= 1) {
            bound = candidate;
            candidate = left_side ? previous_part(bound) : next_part(bound);
        }
    }
};

bool adjacent_by_scan(const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v) {
    if (u == NetworKit::none || v == NetworKit::none) {
        return false;
    }
    if (graph.degree(u) > graph.degree(v)) {
        std::swap(u, v);
    }
    for (auto w : graph.neighborRange(u)) {
        if (w == v) {
            return true;
        }
    }
    return false;
}

void split_origin_part(PartitionRefinement &P, const NetworKit::Graph &graph) {
    auto H = P.vertex_part[P.origin];
    if (P.parts[H].vertices.size() <= 1) {
        return;
    }
    P.left = P.right = H;
    if (P.parts[H].pivot == NetworKit::none) {
        P.parts[H].pivot = P.parts[H].vertices.front();
    }
    auto neighbours = P.add_part_after(H);
    for (auto v : graph.neighborRange(P.origin)) {
        if (P.vertex_part[v] == H) {
            P.move_vertex(v, neighbours);
        }
    }
    if (P.parts[neighbours].vertices.empty()) {
        P.erase_part(neighbours);
    } else {
        P.active_parts++;
        P.unused_parts.push_back(neighbours);
    }

    if (P.parts[H].vertices.size() > 1) {
        auto singleton = P.add_part_after(H);
        P.move_vertex(P.origin, singleton);
        P.parts[singleton].pivot = P.origin;
        P.active_parts++;
        P.unused_parts.push_back(H);
        P.left = P.right = singleton;
    }
    P.check_bound(P.left, true);
    P.check_bound(P.right, false);
}

void refine_with_pivot(
        PartitionRefinement &P, const NetworKit::Graph &graph, NetworKit::node H) {
    if (!P.parts[H].active || P.parts[H].vertices.empty()) {
        return;
    }
    if (P.parts[H].pivot == NetworKit::none) {
        P.parts[H].pivot = P.parts[H].vertices.front();
    }
    auto pivot = P.parts[H].pivot;
    for (auto v : graph.neighborRange(pivot)) {
        auto part = P.vertex_part[v];
        if (part != H && part != NetworKit::none && P.parts[part].active) {
            P.parts[part].amount++;
        }
    }
    for (auto v : graph.neighborRange(pivot)) {
        auto part = P.vertex_part[v];
        if (part == H || part == NetworKit::none || !P.parts[part].active
                || P.parts[part].vertices.size() <= P.parts[part].amount) {
            continue;
        }
        if (P.parts[part].division == NetworKit::none) {
            if (P.parts[part].vertices.size() == 1) {
                continue;
            }
            P.parts[part].division = P.add_part_after(part);
            P.active_parts++;
        }
        auto division = P.parts[part].division;
        if (P.parts[part].pivot == v) {
            P.unused_parts.push_back(part);
            P.parts[division].pivot = v;
        }
        P.move_vertex(v, division);
        P.parts[part].amount--;
    }

    for (auto v : graph.neighborRange(pivot)) {
        auto part = P.vertex_part[v];
        if (part == NetworKit::none || !P.parts[part].active) {
            continue;
        }
        P.parts[part].amount = 0;
        auto previous = P.previous_part(part);
        if (part != H && part != P.order.front() && previous != NetworKit::none
                && P.parts[previous].division == part) {
            P.parts[previous].division = NetworKit::none;
            P.parts[previous].amount = 0;
            if (P.parts[previous].vertices.empty()) {
                P.erase_part(previous);
                P.active_parts--;
            }
            if (P.parts[part].pivot == NetworKit::none) {
                P.unused_parts.push_back(part);
            }
        }
    }
}

void choose_new_origin(PartitionRefinement &P, const NetworKit::Graph &graph) {
    auto origin_part = P.vertex_part[P.origin];
    auto left_pivot = P.parts[P.left].pivot;
    auto right_pivot = P.parts[P.right].pivot;
    P.origin = P.right == origin_part || (P.left != origin_part
            && adjacent_by_scan(graph, left_pivot, right_pivot)) ? left_pivot : right_pivot;
}

std::vector<NetworKit::node> factorizing_permutation(const NetworKit::Graph &graph) {
    PartitionRefinement P(graph);
    while (P.active_parts < graph.numberOfNodes()) {
        split_origin_part(P, graph);
        while (!P.unused_parts.empty()) {
            auto H = P.unused_parts.front();
            P.unused_parts.pop_front();
            refine_with_pivot(P, graph, H);
            P.check_bound(P.left, true);
            P.check_bound(P.right, false);
        }
        choose_new_origin(P, graph);
    }

    std::vector<NetworKit::node> permutation;
    permutation.reserve(graph.numberOfNodes());
    for (auto part : P.order) {
        if (!P.parts[part].vertices.empty()) {
            permutation.push_back(P.parts[part].vertices.front());
        }
    }
    return permutation;
}

bool same_neighborhood(
        const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v,
        std::vector<NetworKit::count> &used, NetworKit::count marker, bool adjacent) {
    if (u == NetworKit::none || v == NetworKit::none || graph.degree(u) != graph.degree(v)) {
        return false;
    }
    NetworKit::count common = 0;
    bool saw_v = false, saw_u = false;
    for (auto w : graph.neighborRange(u)) {
        if (adjacent && w == v) {
            saw_v = true;
            continue;
        }
        used[w] = marker;
        common++;
    }
    for (auto w : graph.neighborRange(v)) {
        if (adjacent && w == u) {
            saw_u = true;
            continue;
        }
        if (used[w] != marker) {
            return false;
        }
        common--;
    }
    return common == 0 && (!adjacent || (saw_v && saw_u));
}

NodeType twins(
        const NetworKit::Graph &graph, NetworKit::node u, NetworKit::node v,
        std::vector<NetworKit::count> &used, NetworKit::count &marker) {
    if (same_neighborhood(graph, u, v, used, ++marker, true)) {
        return NodeType::COMPLEMENT_NODE;
    }
    if (same_neighborhood(graph, u, v, used, ++marker, false)) {
        return NodeType::UNION_NODE;
    }
    return NodeType::UNKNOWN;
}

void build_cotree(
        Cotree &cotree, std::vector<TwinOrder> order,
        NetworKit::count upper_node_id_bound) {
    std::reverse(order.begin(), order.end());
    cotree.clear();
    cotree.prepared = true;
    if (order.empty()) {
        return;
    }

    const NetworKit::node n = order.size() + upper_node_id_bound - 1;
    cotree.reserve(2 * n);
    auto root = cotree.add(NodeType::UNION_NODE);
    cotree.setRoot(root);

    std::vector<NetworKit::node> leaf(upper_node_id_bound, NetworKit::none);
    auto first = std::get<0>(order[0]);
    leaf[first] = cotree.add(NodeType::LEAF, first);
    cotree.addChild(root, leaf[first]);

    for (NetworKit::index i = 1; i < order.size(); ++i) {
        auto [removed, existing, node_type] = order[i];
        auto parent = cotree.getNode(leaf[existing]).parent;
        auto internal = cotree.add(node_type);

        cotree.replaceChild(parent, leaf[existing], internal);
        leaf[removed] = cotree.add(NodeType::LEAF, removed);
        cotree.addChild(internal, leaf[existing]);
        cotree.addChild(internal, leaf[removed]);
    }
}

}  // namespace

void HabibPaulCographRecognition::run() {
    hasRun = true;
    cotree.clear();
    if (graph.numberOfNodes() == 0) {
        is_cograph = State::COGRAPH;
        return;
    }

    auto permutation = factorizing_permutation(graph);
    NetworKit::Graph working_graph(graph);
    std::list<NetworKit::node> remaining;
    remaining.push_back(NetworKit::none);
    remaining.insert(remaining.end(), permutation.begin(), permutation.end());
    remaining.push_back(NetworKit::none);

    std::vector<TwinOrder> order;
    std::vector<NetworKit::count> used(graph.upperNodeIdBound(), 0);
    NetworKit::count marker = 0, removed = 0;
    auto z = std::next(remaining.begin());
    while (z != std::prev(remaining.end())) {
        auto previous = std::prev(z);
        auto next = std::next(z);
        auto twin_type = twins(working_graph, *z, *previous, used, marker);
        if (twin_type != NodeType::UNKNOWN) {
            order.emplace_back(*previous, *z, twin_type);
            working_graph.removeNode(*previous);
            remaining.erase(previous);
            removed++;
            continue;
        }
        twin_type = twins(working_graph, *z, *next, used, marker);
        if (twin_type != NodeType::UNKNOWN) {
            order.emplace_back(*z, *next, twin_type);
            working_graph.removeNode(*z);
            z = remaining.erase(z);
            removed++;
        } else {
            ++z;
        }
    }

    if (removed == graph.numberOfNodes() - 1) {
        is_cograph = State::COGRAPH;
        order.emplace_back(*std::next(remaining.begin()), NetworKit::none, NodeType::UNKNOWN);
        build_cotree(cotree, order, graph.upperNodeIdBound());
    } else {
        is_cograph = State::NOT_COGRAPH;
    }
}

}  // namespace Koala
