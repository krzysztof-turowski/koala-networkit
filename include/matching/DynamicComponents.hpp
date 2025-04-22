#pragma once

#include <networkit/graph/Graph.hpp>
#include <Eigen/Core>


namespace Koala {
    class DynamicComponents {
        public:
        DynamicComponents();
        DynamicComponents(const NetworKit::Graph& graph, const Eigen::MatrixXd& AG);
        void addVertex();
        void addVertex(int count);
        void removeVertex(int idx);
        void addEdge(int u, int v);
        void removeEdge(int u, int v);
        bool isConnected(int u, int v) const;
        int getComponentSize(int u) const;
        int size() const;
        int getLabel(int i) const;

        NetworKit::Graph G;
        Eigen::MatrixXd AG;
    };
}
