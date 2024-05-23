#pragma once

#include "VertexColoring.hpp"
#include "recognition/CographAlg.hpp"

namespace Koala {

    class CographVertexColoring : public VertexColoring {
    private:
        void end_of_coloring(NetworKit::count v);

        void subtree_colors(NetworKit::count v);

        std::vector<NetworKit::count> color, number_of_colors;

        Koala::CographRecognition *recognition = new Koala::CographRecognition(graph.value());
    public:
        using VertexColoring::VertexColoring;

        void run();

        bool checkColoring();

        const std::map<NetworKit::node, int> &getColoring() const;
    };

}
