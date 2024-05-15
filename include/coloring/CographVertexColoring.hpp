#pragma once

#include "VertexColoring.hpp"
#include "recognition/CographAlg.hpp"

namespace Koala {

    class CographVertexColoring : public VertexColoring {
    public:
        using VertexColoring::VertexColoring;

        std::vector<NetworKit::count> color, number_of_colors;

        Koala::CographRecognition *recognition = new Koala::CographRecognition(graph.value());

        void run();

        void EndOfColoring(NetworKit::count v);

        void SubtreeColors(NetworKit::count v);

        NetworKit::count GetColor(NetworKit::count i);

        bool CheckColoring();
    };

}
