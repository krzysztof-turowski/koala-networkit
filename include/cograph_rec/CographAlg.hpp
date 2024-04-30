#pragma once

#include <optional>
#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#include "cograph_rec/Cotree.hpp"
#include "cograph_rec/FactorizingPermutation.hpp"
#include "cograph_rec/Part.hpp"

namespace Koala {
    class CographRecognition : public NetworKit::Algorithm {
    public:
        enum class State {
            UNKNOWN,
            COGRAPH,
            IS_NOT_COGRAPH
        };
        State status;
        NetworKit::Graph *graph, *original_graph;
        Cotree *cotree;
        FactorizingPermutation *permutation;
        long long num_of_parts, num_of_nodes;

        std::vector<std::pair<std::pair<long long, long long>, long long> > order;

        std::list<part *> unused_parts;

        bool IsCograph() const;

        CographRecognition::State GetState() const;

        explicit CographRecognition(NetworKit::Graph &graph);

        ~CographRecognition();

        void run();

        void Clear();
    };
} /* namespace Koala */
