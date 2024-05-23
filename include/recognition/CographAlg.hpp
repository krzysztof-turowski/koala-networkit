#pragma once

#include <optional>
#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#include "structures/Cotree.hpp"
#include "Cograph/FactorizingPermutation.hpp"
#include "Cograph/Part.hpp"
#include "Cograph/Twins.hpp"

namespace Koala {
    class CographRecognition : public NetworKit::Algorithm {
    public:
        enum class State {
            UNKNOWN,
            COGRAPH,
            NOT_COGRAPH
        };

        bool isCograph() const;

        CographRecognition::State getState() const;

        CographRecognition(NetworKit::Graph &graph);

        ~CographRecognition();

        void run();

        Cotree cotree;

        private:

        void Clear();

        State status;
        NetworKit::Graph graph, original_graph;
        Twins T;

        FactorizingPermutation permutation;
        NetworKit::count num_of_parts, num_of_nodes;

        std::vector<std::pair<std::pair<NetworKit::count, NetworKit::count>, NetworKit::count> > order;
        part *H;
        std::list<part *> unused_parts;
    };
}
