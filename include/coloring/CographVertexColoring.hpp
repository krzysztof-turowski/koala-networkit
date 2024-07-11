#pragma once

#include "VertexColoring.hpp"
#include "recognition/CographRecognition.hpp"

namespace Koala {

    class CographVertexColoring : public VertexColoring {
    private:
        void end_of_coloring();

        void subtree_colors();

        std::vector<NetworKit::count> color, number_of_colors;

        Koala::Cotree &cotree;

        std::stack<int> st;

        std::vector<bool> used;
    public:
        CographVertexColoring(NetworKit::Graph Graph, Koala::Cotree &Cotree) : VertexColoring(Graph), cotree(Cotree) {
        }

        void run();

        bool checkColoring();
    };
} /* namespace Koala */

