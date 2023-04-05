#pragma once

#include "MinimalSpanningTree.hpp"
#include "utils.hpp"


namespace MST {
    class KruskalMST final : public MinimalSpanningTree {
    public:
        explicit KruskalMST(const Graph& G, bool storeOriginalGraph = false);

        ~KruskalMST() override = default;
    };
}
