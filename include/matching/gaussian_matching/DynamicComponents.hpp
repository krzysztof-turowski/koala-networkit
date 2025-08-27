#pragma once

#include <networkit/graph/Graph.hpp>
#include <eigen3/Eigen/Core>

namespace Koala {
    class DynamicComponents {
        public:
        DynamicComponents(const NetworKit::Graph& G);

        void addEdge(int u, int v);
        void removeEdge(int u, int v);
        bool isConnected(int u, int v) const;
        int getComponentSize(int v) const;

        NetworKit::Graph G;
    };
}
