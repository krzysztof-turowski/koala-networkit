#pragma once

#include <flow/MinimumCostFlow.hpp>

namespace Koala{
    class SuccesiveApproximation final : public MinimumCostFlow {
        using MinimumCostFlow::MinimumCostFlow;
    };
}