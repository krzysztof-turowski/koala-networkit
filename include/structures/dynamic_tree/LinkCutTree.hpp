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

#include <networkit/graph/Graph.hpp>

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

} /* namespace Koala */
