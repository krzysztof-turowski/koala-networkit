#include "structures/Cotree.hpp"

namespace Koala {

Cotree::Cotree() {
    prepared = false;
    root = NetworKit::none;
}

NetworKit::node Cotree::getRoot() const {
    return root;
}

NetworKit::node Cotree::upperNodeIdBound() const {
    return nodes.size();
}

void Cotree::reserve(NetworKit::count n) {
    nodes.reserve(n);
}

NetworKit::node Cotree::add(NodeType type, NetworKit::node number) {
    NetworKit::node node = nodes.size();
    nodes.emplace_back(NetworKit::none, NetworKit::none, NetworKit::none);
    nodes[node].type = type;
    if (type == NodeType::LEAF) {
        nodes[node].number = number;
    }
    return node;
}

void Cotree::clear() {
    nodes.clear();
    root = NetworKit::none;
    prepared = false;
}

void Cotree::setRoot(NetworKit::node node) {
    root = node;
    if (node != NetworKit::none) {
        nodes[node].parent = NetworKit::none;
        nodes[node].previous_sibling = NetworKit::none;
        nodes[node].next_sibling = NetworKit::none;
    }
}

void Cotree::addChild(NetworKit::node parent, NetworKit::node child) {
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
    nodes[parent].size++;
}

void Cotree::removeChild(NetworKit::node parent, NetworKit::node child) {
    if (nodes[child].parent != parent) {
        return;
    }
    NetworKit::node previous = nodes[child].previous_sibling;
    NetworKit::node next = nodes[child].next_sibling;
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
    nodes[parent].size--;
}

void Cotree::replaceChild(
        NetworKit::node parent, NetworKit::node old_child, NetworKit::node new_child) {
    NetworKit::node previous = nodes[old_child].previous_sibling;
    NetworKit::node next = nodes[old_child].next_sibling;
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

void Cotree::moveChildToFront(NetworKit::node parent, NetworKit::node child) {
    if (nodes[child].parent != parent || nodes[parent].first_child == child) {
        return;
    }
    removeChild(parent, child);
    addChild(parent, child);
}

} /* namespace Koala */
