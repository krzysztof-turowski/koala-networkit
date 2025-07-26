#pragma once

#include <networkit/graph/Graph.hpp>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <matching/gaussian_matching/utils.hpp>
#include <set>

typedef std::set<std::pair<int, int>> Matching;

namespace Koala {
    class BipartiteGaussianMatching {
        friend class BpTest;

    public:
        BipartiteGaussianMatching(const NetworKit::Graph& G);
        void run();
        Matching getMatching();

        // private:
        NetworKit::Graph G;
        MatZp AG;
        Matching M;

        std::vector<int> U, V;
        std::vector<int> bpIdx;
        std::vector<int> oldIdx;
    };
}
