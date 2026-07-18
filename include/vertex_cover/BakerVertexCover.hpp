#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

#include "techniques/BakerKOuterplanarGraphScheme.hpp"

namespace Koala {

/**
 * Problem adapter for minimum vertex cover in the Baker scheme.
 *
 * State 0 excludes a vertex and state 1 includes it. Every original edge must
 * have at least one included endpoint. Table values are minimized.
 */
class BakerVertexCover {
 public:
    using Value = std::int64_t;
    using Solution = std::set<NetworKit::node>;

    static constexpr std::size_t EXCLUDED = 0;
    static constexpr std::size_t INCLUDED = 1;

    bool isMaximization() const {
        return false;
    }

    std::size_t stateCount() const {
        return 2;
    }

    Value infeasibleValue() const {
        return std::numeric_limits<Value>::max() / 4;
    }

    Value identityValue() const {
        return 0;
    }

    bool isInfeasible(Value value) const {
        return value == infeasibleValue();
    }

    Value vertexValue(NetworKit::node, std::size_t state) const {
        return state == INCLUDED ? 1 : 0;
    }

    bool isValidEdge(
            NetworKit::node, std::size_t first_state,
            NetworKit::node, std::size_t second_state) const {
        return first_state == INCLUDED || second_state == INCLUDED;
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

    std::optional<Value> directValue(
            const NetworKit::Graph &graph,
            const BakerSliceBoundary &boundary,
            const std::vector<std::size_t> &states) const {
        if (states.size() != boundary.size()) {
            return std::nullopt;
        }
        std::vector<NetworKit::node> vertices = boundary.left;
        vertices.insert(
            vertices.end(), boundary.right.begin(), boundary.right.end());

        Value value = identityValue();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            bool duplicate = false;
            for (std::size_t j = 0; j < i; ++j) {
                if (vertices[i] == vertices[j]) {
                    duplicate = true;
                    if (states[i] != states[j]) {
                        return std::nullopt;
                    }
                }
            }
            if (graph.hasEdge(vertices[i], vertices[i])
                && !isValidEdge(
                    vertices[i], states[i], vertices[i], states[i])) {
                return std::nullopt;
            }
            if (!duplicate) {
                value = combineValues(
                    value, vertexValue(vertices[i], states[i]));
            }
        }
        for (std::size_t first = 0; first < vertices.size(); ++first) {
            for (std::size_t second = first + 1;
                 second < vertices.size(); ++second) {
                if (vertices[first] != vertices[second]
                    && graph.hasEdge(vertices[first], vertices[second])
                    && !isValidEdge(
                        vertices[first], states[first],
                        vertices[second], states[second])) {
                    return std::nullopt;
                }
            }
        }
        return value;
    }

    std::size_t adjustmentCount() const {
        return 1;
    }

    bool adjustInputStates(
            const NetworKit::Graph &graph,
            const BakerSliceBoundary &boundary,
            NetworKit::node x, NetworKit::node y,
            const std::vector<std::size_t> &output_states,
            std::size_t transition,
            std::vector<std::size_t> &input_states) const {
        const std::size_t level = boundary.left.size();
        if (transition != 0 || level == 0
            || boundary.right.size() != level
            || output_states.size() != 2 * level) {
            return false;
        }
        input_states = output_states;
        if (x == y) {
            return output_states[0] == output_states[level]
                && (!graph.hasEdge(x, x)
                    || isValidEdge(
                        x, output_states[0], x, output_states[0]));
        }
        return !graph.hasEdge(x, y)
            || isValidEdge(
                x, output_states[0], y, output_states[level]);
    }

    Value adjustValue(
            Value input_value, NetworKit::node x, NetworKit::node y,
            const std::vector<std::size_t> &output_states) const {
        if (x == y) {
            return removeDuplicate(input_value, x, output_states[0]);
        }
        return input_value;
    }

    bool splitMergeState(
            std::size_t merge_state, std::size_t &first_state,
            std::size_t &second_state) const {
        if (merge_state >= stateCount()) {
            return false;
        }
        first_state = merge_state;
        second_state = merge_state;
        return true;
    }

    bool mergeInputStates(
            NetworKit::node, NetworKit::node, NetworKit::node,
            std::size_t output_left_state, std::size_t merge_state,
            std::size_t output_right_state,
            std::size_t &first_left_state,
            std::size_t &first_middle_state,
            std::size_t &second_middle_state,
            std::size_t &second_right_state) const {
        if (output_left_state >= stateCount()
            || output_right_state >= stateCount()) {
            return false;
        }
        first_left_state = output_left_state;
        second_right_state = output_right_state;
        return splitMergeState(
            merge_state, first_middle_state, second_middle_state);
    }

    bool canonicalMergeOutputStates(
            NetworKit::node, NetworKit::node,
            std::size_t left_state, std::size_t right_state,
            std::size_t &output_left_state,
            std::size_t &output_right_state) const {
        output_left_state = left_state;
        output_right_state = right_state;
        return left_state < stateCount() && right_state < stateCount();
    }

    Value mergeValue(
            Value first, Value second,
            const std::vector<NetworKit::node> &middle_boundary,
            const std::vector<std::size_t> &middle_states) const {
        if (middle_boundary.size() != middle_states.size()) {
            return infeasibleValue();
        }
        Value value = combineValues(first, second);
        for (std::size_t i = 0; i < middle_boundary.size(); ++i) {
            value = removeDuplicate(
                value, middle_boundary[i], middle_states[i]);
        }
        return value;
    }

    bool canContract(NetworKit::node, std::size_t state) const {
        return state < stateCount();
    }

    bool extendInputStates(
            const NetworKit::Graph &graph, NetworKit::node vertex,
            const BakerSliceBoundary &input_boundary,
            const std::vector<std::size_t> &output_states,
            std::vector<std::size_t> &input_states) const {
        const std::size_t input_level = input_boundary.left.size();
        const std::size_t output_level = input_level + 1;
        if (input_boundary.right.size() != input_level
            || output_states.size() != 2 * output_level
            || output_states[0] != output_states[output_level]) {
            return false;
        }
        const std::size_t state = output_states[0];
        if (graph.hasEdge(vertex, vertex)
            && !isValidEdge(vertex, state, vertex, state)) {
            return false;
        }

        input_states.resize(2 * input_level);
        for (std::size_t i = 0; i < input_level; ++i) {
            input_states[i] = output_states[i + 1];
            input_states[input_level + i] =
                output_states[output_level + i + 1];
            if (graph.hasEdge(vertex, input_boundary.left[i])
                && !isValidEdge(
                    vertex, state, input_boundary.left[i], input_states[i])) {
                return false;
            }
            if (graph.hasEdge(vertex, input_boundary.right[i])
                && !isValidEdge(
                    vertex, state, input_boundary.right[i],
                    input_states[input_level + i])) {
                return false;
            }
        }
        return true;
    }

    Value extendValue(
            Value input_value, NetworKit::node vertex,
            std::size_t state) const {
        return combineValues(input_value, vertexValue(vertex, state));
    }

    bool isFinalAssignment(
            const BakerSliceBoundary &boundary,
            const std::vector<std::size_t> &states) const {
        return states.size() == boundary.size()
            && std::all_of(
                states.begin(), states.end(), [this](std::size_t state) {
                    return state < stateCount();
                });
    }

    bool better(Value first, Value second) const {
        return first < second;
    }

    void appendToSolution(
            NetworKit::node vertex, std::size_t state,
            Solution &solution) const {
        if (state == INCLUDED) {
            solution.insert(vertex);
        }
    }

    void finalizeSolution(Solution &) const { }
};

}  // namespace Koala
