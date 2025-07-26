#pragma once

#include <networkit/graph/Graph.hpp>
#include <set>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <matching/utils.hpp>

typedef std::set<std::pair<int, int>> Matching;

namespace Koala {
    class GeneralGaussianMatching {
    public:
        GeneralGaussianMatching(const NetworKit::Graph& G);
        void run();
        Matching getMatching();

        // private:
        NetworKit::Graph G;
        MatZp AG;
        Matching M;

        std::vector<int> oldIdx;
    };
}
