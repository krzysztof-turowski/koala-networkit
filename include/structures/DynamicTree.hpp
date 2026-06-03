#pragma once

#include <networkit/graph/Graph.hpp>

namespace Koala {

template <typename Value>
class LinkCutTree;

template <typename Value>
class DynamicTree {
 public:
    virtual ~DynamicTree() = default;

    virtual void link(NetworKit::node, NetworKit::node, Value) = 0;
    virtual void cut(NetworKit::node, NetworKit::node) = 0;
    virtual NetworKit::node findRoot(NetworKit::node) = 0;
};

template <class T>
using DefaultDynamicTree = LinkCutTree<T>;

} /* namespace Koala */
