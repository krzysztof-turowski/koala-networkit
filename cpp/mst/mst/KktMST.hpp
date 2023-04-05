#pragma once

#include "MinimalSpanningTree.hpp"
#include "MSTV.hpp"

namespace MST {
    class KktMST final : public MinimalSpanningTree {
    public:
        explicit KktMST(const Graph& G, bool storeOriginalGraph = false);

        ~KktMST() override = default;
    };
}
