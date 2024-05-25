#pragma once

#include "recognition/CographAlg.hpp"
#include "Pathwidth.hpp"

namespace Koala {
    class CographPathwidth : public Pathwidth{
    private:
        NetworKit::count pathwidth(NetworKit::count v);

        NetworKit::count subtree_size(NetworKit::count v);
    public:

        Koala::Cotree &cotree;

        CographPathwidth(NetworKit::Graph &Graph,Koala::Cotree &CoTree) : Pathwidth(Graph),cotree(CoTree){

        }

        void run();
    };
}
