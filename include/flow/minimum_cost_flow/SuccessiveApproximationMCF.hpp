#pragma once

#include <flow/MinimumCostFlow.hpp>

namespace Koala{
    class SuccessiveApproximationMCF final : public MinimumCostFlow {
        using MinimumCostFlow::MinimumCostFlow;
    };
}