#pragma once

#include <graph/GraphTools.hpp>

namespace Koala {
    class CoNode {
    public:
        NetworKit::count left_son, right_son, parent, type, size;

        CoNode(NetworKit::count l, NetworKit::count r, NetworKit::count p) {
            left_son = l;
            right_son = r;
            parent = p;
            type = 0;
            size = 0;
        }
    };

    class Cotree {
    public:
        NetworKit::Graph *graph;
        std::vector<CoNode> nodes;
        std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > order;

        explicit Cotree(NetworKit::Graph &Graph);

        bool prepared;

        void BuildTree();

        void SetOrder(std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > a) {
            order = a;
        }
    };
}