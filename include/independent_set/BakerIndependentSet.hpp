#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * Problem adapter for maximum independent set in the Baker scheme.
 *
 * State 0 excludes a vertex and state 1 includes it. Table values are
 * maximized, and an original edge forbids assigning state 1 to both endpoints.
 */
class BakerIndependentSet {
 public:
    using Value = std::int64_t;
    using Solution = std::vector<NetworKit::node>;

    std::size_t stateCount() const {
        return 2;
    }

    Value infeasibleValue() const {
        return std::numeric_limits<Value>::lowest() / 4;
    }

    Value identityValue() const {
        return 0;
    }

    bool isInfeasible(Value value) const {
        return value == infeasibleValue();
    }

    Value vertexValue(NetworKit::node, std::size_t state) const {
        return state == 1 ? 1 : 0;
    }

    bool isValidEdge(
        NetworKit::node, std::size_t first_state,
        NetworKit::node, std::size_t second_state) const {
        return first_state == 0 || second_state == 0;
    }

    Value combineValues(Value first, Value second) const {
        if (isInfeasible(first) || isInfeasible(second)) {
            return infeasibleValue();
        }
        return first + second;
    }

    Value removeDuplicate(
            Value value, NetworKit::node vertex, std::size_t state) const {
        return value - vertexValue(vertex, state);
    }

    bool better(Value first, Value second) const {
        return first > second;
    }

    void appendToSolution(
        NetworKit::node vertex, std::size_t state, Solution &solution) const {
        if (state == 1) {
            solution.push_back(vertex);
        }
    }

    void finalizeSolution(Solution &solution) const {
        std::sort(solution.begin(), solution.end());
        solution.erase(std::unique(solution.begin(), solution.end()), solution.end());
    }
};

}  // namespace Koala
