/*
 * DynamicTree.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <memory>
#include <optional>
#include <unordered_set>
#include <vector>

#include <structures/DynamicTree.hpp>

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

template <typename Value>
class LinkCutTree final : public DynamicTree<Value> {
 public:
    LinkCutTree() = default;

    void initialize(NetworKit::count);

    NetworKit::node findParent(NetworKit::node);
    NetworKit::node findRoot(NetworKit::node) override;
    const std::unordered_set<NetworKit::node>& findChildren(NetworKit::node) const;

    std::optional<Value> getValue(NetworKit::node);
    void addValue(NetworKit::node, Value);

    void link(NetworKit::node, NetworKit::node, Value) override;
    void cut(NetworKit::node, NetworKit::node) override;

    NetworKit::Edge findSaturatedEdge(NetworKit::node);
    std::optional<Value> getMinimumPathResidualCapacity(NetworKit::node);

 private:
    struct Node {
        explicit Node(NetworKit::node id) : id(id) {}

        NetworKit::node id;
        Node *left = nullptr, *right = nullptr, *parent = nullptr;
        std::optional<Value> value, minimum;
        Value lazy = 0;
    };

    std::vector<std::unique_ptr<Node>> nodes;
    std::vector<std::unordered_set<NetworKit::node>> children;

    static bool isAuxiliaryRoot(const Node*);
    static void add(Node*, Value);
    static void push(Node*);
    static void pull(Node*);
    static void rotate(Node*);
    static void splay(Node*);
    static void access(Node*);
    static Node* findRoot(Node*);
    static Node* findParent(Node*);
    static Node* findRightmostAtMost(Node*, Value);

    Node* get(NetworKit::node);
};

template <typename Value>
bool LinkCutTree<Value>::isAuxiliaryRoot(const Node *node) {
    return node->parent == nullptr
        || (node->parent->left != node && node->parent->right != node);
}

template <typename Value>
void LinkCutTree<Value>::add(Node *node, Value delta) {
    if (node == nullptr || !node->minimum.has_value()) {
        return;
    }
    if (node->value.has_value()) {
        *node->value += delta;
    }
    *node->minimum += delta;
    node->lazy += delta;
}

template <typename Value>
void LinkCutTree<Value>::push(Node *node) {
    if (node->lazy == 0) {
        return;
    }
    add(node->left, node->lazy);
    add(node->right, node->lazy);
    node->lazy = 0;
}

template <typename Value>
void LinkCutTree<Value>::pull(Node *node) {
    node->minimum = node->value;
    for (const auto *child : {node->left, node->right}) {
        if (child != nullptr && child->minimum.has_value()
                && (!node->minimum.has_value() || *child->minimum < *node->minimum)) {
            node->minimum = child->minimum;
        }
    }
}

template <typename Value>
void LinkCutTree<Value>::rotate(Node *node) {
    Node *parent = node->parent;
    Node *grandparent = parent->parent;
    if (!isAuxiliaryRoot(parent)) {
        if (grandparent->left == parent) {
            grandparent->left = node;
        } else {
            grandparent->right = node;
        }
    }
    node->parent = grandparent;

    if (parent->left == node) {
        parent->left = node->right;
        if (node->right != nullptr) {
            node->right->parent = parent;
        }
        node->right = parent;
    } else {
        parent->right = node->left;
        if (node->left != nullptr) {
            node->left->parent = parent;
        }
        node->left = parent;
    }
    parent->parent = node;
    pull(parent);
    pull(node);
}

template <typename Value>
void LinkCutTree<Value>::splay(Node *node) {
    std::vector<Node *> ancestors{node};
    for (Node *parent = node; !isAuxiliaryRoot(parent); parent = parent->parent) {
        ancestors.push_back(parent->parent);
    }
    for (auto it = ancestors.rbegin(); it != ancestors.rend(); ++it) {
        push(*it);
    }

    while (!isAuxiliaryRoot(node)) {
        Node *parent = node->parent;
        Node *grandparent = parent->parent;
        if (!isAuxiliaryRoot(parent)) {
            if ((parent->left == node) == (grandparent->left == parent)) {
                rotate(parent);
            } else {
                rotate(node);
            }
        }
        rotate(node);
    }
}

template <typename Value>
void LinkCutTree<Value>::access(Node *node) {
    Node *previous = nullptr;
    for (Node *current = node; current != nullptr; current = current->parent) {
        splay(current);
        current->right = previous;
        if (previous != nullptr) {
            previous->parent = current;
        }
        pull(current);
        previous = current;
    }
    splay(node);
}

template <typename Value>
typename LinkCutTree<Value>::Node* LinkCutTree<Value>::findRoot(Node *node) {
    access(node);
    while (node->left != nullptr) {
        push(node);
        node = node->left;
    }
    splay(node);
    return node;
}

template <typename Value>
typename LinkCutTree<Value>::Node* LinkCutTree<Value>::findParent(Node *node) {
    access(node);
    if (node->left == nullptr) {
        return nullptr;
    }
    node = node->left;
    push(node);
    while (node->right != nullptr) {
        node = node->right;
        push(node);
    }
    splay(node);
    return node;
}

template <typename Value>
typename LinkCutTree<Value>::Node* LinkCutTree<Value>::findRightmostAtMost(
    Node *node, Value limit) {
    push(node);
    if (node->right != nullptr && node->right->minimum.has_value()
            && *node->right->minimum <= limit) {
        return findRightmostAtMost(node->right, limit);
    }
    if (node->value.has_value() && *node->value <= limit) {
        return node;
    }
    if (node->left != nullptr && node->left->minimum.has_value()
            && *node->left->minimum <= limit) {
        return findRightmostAtMost(node->left, limit);
    }
    return nullptr;
}

template <typename Value>
typename LinkCutTree<Value>::Node* LinkCutTree<Value>::get(NetworKit::node v) {
    return nodes[v].get();
}

template <typename Value>
void LinkCutTree<Value>::initialize(NetworKit::count n) {
    nodes.clear();
    nodes.reserve(n);
    for (NetworKit::node v = 0; v < n; ++v) {
        nodes.push_back(std::make_unique<Node>(v));
    }
    children.assign(n, {});
}

template <typename Value>
NetworKit::node LinkCutTree<Value>::findParent(NetworKit::node v) {
    const auto *parent = findParent(get(v));
    return parent == nullptr ? NetworKit::none : parent->id;
}

template <typename Value>
NetworKit::node LinkCutTree<Value>::findRoot(NetworKit::node v) {
    return findRoot(get(v))->id;
}

template <typename Value>
const std::unordered_set<NetworKit::node>& LinkCutTree<Value>::findChildren(
    NetworKit::node v) const {
    return children[v];
}

template <typename Value>
std::optional<Value> LinkCutTree<Value>::getValue(NetworKit::node v) {
    auto *node = get(v);
    access(node);
    return node->value;
}

template <typename Value>
void LinkCutTree<Value>::addValue(NetworKit::node v, Value delta) {
    auto *node = get(v);
    access(node);
    add(node, delta);
}

template <typename Value>
void LinkCutTree<Value>::link(NetworKit::node u, NetworKit::node v, Value capacity) {
    auto *child = get(u);
    auto *parent = get(v);
    if (findParent(child) != nullptr || findRoot(parent) == child) {
        return;
    }
    access(child);
    child->parent = parent;
    child->value = capacity;
    pull(child);
    children[v].insert(u);
}

template <typename Value>
void LinkCutTree<Value>::cut(NetworKit::node u, NetworKit::node v) {
    auto *child = get(u);
    auto *parent = findParent(child);
    if (parent == nullptr || parent->id != v) {
        return;
    }
    access(child);
    child->left->parent = nullptr;
    child->left = nullptr;
    child->value.reset();
    pull(child);
    children[v].erase(u);
}

template <typename Value>
NetworKit::Edge LinkCutTree<Value>::findSaturatedEdge(NetworKit::node v) {
    auto *node = get(v);
    access(node);
    auto *bottleneck = findRightmostAtMost(node, 0);
    if (bottleneck == nullptr) {
        return NetworKit::Edge();
    }
    auto *parent = findParent(bottleneck);
    return parent == nullptr ? NetworKit::Edge()
                             : NetworKit::Edge(bottleneck->id, parent->id);
}

template <typename Value>
std::optional<Value> LinkCutTree<Value>::getMinimumPathResidualCapacity(
    NetworKit::node v) {
    auto *node = get(v);
    access(node);
    return node->minimum;
}

} /* namespace Koala */
