#pragma once

#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {
    class Pathwidth : public NetworKit::Algorithm {
    public:
        virtual void run() = 0;

        NetworKit::count getPathwidthSize();

        Pathwidth(const NetworKit::Graph &graph);

    protected:
        NetworKit::count width = 0;
        std::optional<NetworKit::Graph> graph;

    };
}
