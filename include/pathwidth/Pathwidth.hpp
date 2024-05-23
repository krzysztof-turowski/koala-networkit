#pragma once

#include "recognition/CographAlg.hpp"

namespace Koala {
    class Pathwidth {
    private:
        NetworKit::count pathwidth(NetworKit::count v);

        NetworKit::count subtree_size(NetworKit::count v);
    public:

        Koala::Cotree &cotree;

        Pathwidth(Koala::Cotree &CoTree) : cotree(CoTree){

        }

        NetworKit::count width = 0;

        void run();

    };
}
