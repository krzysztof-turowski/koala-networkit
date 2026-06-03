#include <graph/GraphTools.hpp>

#include "structures/Cotree.hpp"

namespace Koala {

Cotree::Cotree() {
    graph = nullptr;
    prepared = false;
    root = NetworKit::none;
}

Cotree::Cotree(NetworKit::Graph &Graph) {
    graph = &Graph;
    prepared = false;
    root = NetworKit::none;
}

NetworKit::count Cotree::getRoot() const {
    return root;
}

NetworKit::count Cotree::upperNodeIdBound() const {
    return nodes.size();
}

void Cotree::reserve(NetworKit::count n) {
    nodes.reserve(n);
}

NetworKit::count Cotree::add(NodeType type, int number) {
    NetworKit::count node = nodes.size();
    nodes.emplace_back(NetworKit::none, NetworKit::none, NetworKit::none);
    nodes[node].type = type;
    nodes[node].number = number;
    return node;
}

void Cotree::clear() {
    nodes.clear();
    root = NetworKit::none;
    prepared = false;
}

void Cotree::setRoot(NetworKit::count node) {
    root = node;
    if (node != NetworKit::none) {
        nodes[node].parent = NetworKit::none;
        nodes[node].previous_sibling = NetworKit::none;
        nodes[node].next_sibling = NetworKit::none;
    }
}

void Cotree::addChild(NetworKit::count parent, NetworKit::count child) {
    if (nodes[child].parent != NetworKit::none) {
        removeChild(nodes[child].parent, child);
    }
    nodes[child].parent = parent;
    nodes[child].previous_sibling = NetworKit::none;
    nodes[child].next_sibling = nodes[parent].first_child;
    if (nodes[parent].first_child != NetworKit::none) {
        nodes[nodes[parent].first_child].previous_sibling = child;
    }
    nodes[parent].first_child = child;
    nodes[parent].d++;
}

void Cotree::removeChild(NetworKit::count parent, NetworKit::count child) {
    if (nodes[child].parent != parent) {
        return;
    }
    NetworKit::count previous = nodes[child].previous_sibling;
    NetworKit::count next = nodes[child].next_sibling;
    if (previous != NetworKit::none) {
        nodes[previous].next_sibling = next;
    } else {
        nodes[parent].first_child = next;
    }
    if (next != NetworKit::none) {
        nodes[next].previous_sibling = previous;
    }
    nodes[child].parent = NetworKit::none;
    nodes[child].previous_sibling = NetworKit::none;
    nodes[child].next_sibling = NetworKit::none;
    nodes[parent].d--;
}

void Cotree::replaceChild(
        NetworKit::count parent, NetworKit::count old_child, NetworKit::count new_child) {
    NetworKit::count previous = nodes[old_child].previous_sibling;
    NetworKit::count next = nodes[old_child].next_sibling;
    if (nodes[new_child].parent != NetworKit::none) {
        removeChild(nodes[new_child].parent, new_child);
    }
    nodes[new_child].parent = parent;
    nodes[new_child].previous_sibling = previous;
    nodes[new_child].next_sibling = next;
    if (previous != NetworKit::none) {
        nodes[previous].next_sibling = new_child;
    } else {
        nodes[parent].first_child = new_child;
    }
    if (next != NetworKit::none) {
        nodes[next].previous_sibling = new_child;
    }
    nodes[old_child].parent = NetworKit::none;
    nodes[old_child].previous_sibling = NetworKit::none;
    nodes[old_child].next_sibling = NetworKit::none;
}

void Cotree::moveChildToFront(NetworKit::count parent, NetworKit::count child) {
    if (nodes[child].parent != parent || nodes[parent].first_child == child) {
        return;
    }
    removeChild(parent, child);
    addChild(parent, child);
}

void Cotree::unmarkForNewIteration(NetworKit::count node) {
    nodes[node].marked = Marked::UNMARKED;
    nodes[node].md = 0;
}

void Cotree::mark(NetworKit::count node) {
    nodes[node].marked = Marked::MARKED;
}

void Cotree::unmark(NetworKit::count node) {
    nodes[node].marked = Marked::MARKED_AND_UNMARKED;
}

std::vector<NetworKit::count> Cotree::removeWereMarked(NetworKit::count node) {
    auto child = nodes[node].first_child;
    std::vector<NetworKit::count> removed;
    while (child != NetworKit::none) {
        auto next = nodes[child].next_sibling;
        removed.push_back(child);
        removeChild(node, child);
        child = next;
        if (child == NetworKit::none || nodes[child].marked != Marked::MARKED_AND_UNMARKED) {
            break;
        }
    }
    return removed;
}

void Cotree::removeWereNotMarked(NetworKit::count node) {
    auto child = nodes[node].first_child;
    while (child != NetworKit::none && nodes[child].marked == Marked::MARKED_AND_UNMARKED) {
        child = nodes[child].next_sibling;
    }
    while (child != NetworKit::none) {
        auto next = nodes[child].next_sibling;
        removeChild(node, child);
        child = next;
    }
}

void Cotree::buildTree() {
    reverse(order.begin(), order.end());
    if (order.empty()) {
        nodes.clear();
        root = NetworKit::none;
        prepared = true;
        return;
    }

    NetworKit::count n = order.size() + graph->numberOfNodes() - 1;
    nodes.assign(2 * n, Conode(NetworKit::none, NetworKit::none, NetworKit::none));
    if (n == 0) {
        root = NetworKit::none;
        prepared = true;
        return;
    }

    root = n;
    nodes[root].type = NodeType::UNION_NODE;
    nodes[order[0].first.first].type = NodeType::LEAF;
    nodes[order[0].first.first].parent = root;
    nodes[order[0].first.first].previous_sibling = NetworKit::none;
    nodes[root].first_child = order[0].first.first;
    nodes[root].d = 1;

    for (NetworKit::count i = 1; i < n; i++) {
        const NetworKit::count new_node = n + i;
        const NetworKit::count leaf = order[i].first.first;
        const NetworKit::count existing = order[i].first.second;
        const NetworKit::count parent = nodes[existing].parent;
        const NetworKit::count previous = nodes[existing].previous_sibling;
        const NetworKit::count next = nodes[existing].next_sibling;

        nodes[leaf].type = NodeType::LEAF;
        nodes[leaf].parent = new_node;
        nodes[leaf].previous_sibling = NetworKit::none;

        if (previous == NetworKit::none) {
            nodes[parent].first_child = new_node;
        } else {
            nodes[previous].next_sibling = new_node;
        }

        nodes[new_node].first_child = leaf;
        nodes[new_node].next_sibling = next;
        nodes[new_node].previous_sibling = previous;
        nodes[new_node].parent = parent;
        nodes[new_node].d = 2;
        nodes[new_node].type = order[i].second == 1
            ? NodeType::COMPLEMENT_NODE : NodeType::UNION_NODE;

        if (nodes[new_node].next_sibling != NetworKit::none) {
            nodes[nodes[new_node].next_sibling].previous_sibling = new_node;
        }
        nodes[leaf].next_sibling = existing;
        nodes[existing].parent = new_node;
        nodes[existing].previous_sibling = leaf;
        nodes[existing].next_sibling = NetworKit::none;
    }
    prepared = true;
}

} /* namespace Koala */
