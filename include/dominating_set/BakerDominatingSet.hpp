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
 * Problem adapter for minimum dominating set in the Baker scheme.
 *
 * SELECTED puts the boundary vertex in the dominating set. CERTIFIED records
 * that this slice supplies a selected neighbor. PENDING makes no domination
 * promise and leaves responsibility to another slice. PENDING remains valid
 * even if the represented partial solution happens to dominate the vertex;
 * this relaxation permits a three-way merge split and an O(3^(3k) n) bound.
 * When a shared physical vertex survives on an output boundary, the same
 * three transition codes route no, left, or right certificate ownership.
 * If a cutpoint occurs on both surviving boundaries, a certificate carried
 * by either occurrence certifies the physical vertex on both occurrences.
 */
class BakerDominatingSet {
 public:
    using Value = std::int64_t;
    using Solution = std::set<NetworKit::node>;

    static constexpr std::size_t SELECTED = 0;
    static constexpr std::size_t CERTIFIED = 1;
    static constexpr std::size_t PENDING = 2;

    bool isMaximization() const {
        return false;
    }

    std::size_t stateCount() const {
        return 3;
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
        return state == SELECTED ? 1 : 0;
    }

    bool isValidEdge(
            NetworKit::node, std::size_t,
            NetworKit::node, std::size_t) const {
        return true;
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
            if (states[i] >= stateCount()) {
                return std::nullopt;
            }
            bool duplicate = false;
            for (std::size_t j = 0; j < i; ++j) {
                if (vertices[i] == vertices[j]) {
                    duplicate = true;
                    if (states[i] != states[j]) {
                        return std::nullopt;
                    }
                }
            }
            if (!duplicate) {
                value = combineValues(
                    value, vertexValue(vertices[i], states[i]));
            }
        }
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (states[i] == CERTIFIED
                && !hasSelectedNeighbor(
                    graph, vertices[i], vertices, states)) {
                return std::nullopt;
            }
        }
        return value;
    }

    std::size_t adjustmentCount() const {
        return 3;
    }

    bool adjustInputStates(
            const NetworKit::Graph &graph,
            const BakerSliceBoundary &boundary,
            NetworKit::node x, NetworKit::node y,
            const std::vector<std::size_t> &output_states,
            std::size_t transition,
            std::vector<std::size_t> &input_states) const {
        const std::size_t level = boundary.left.size();
        if (level == 0 || boundary.right.size() != level
            || output_states.size() != 2 * level
            || !validStates(output_states)) {
            return false;
        }
        input_states = output_states;
        if (x != y) {
            if (transition != 0) {
                return false;
            }
            if (graph.hasEdge(x, y)) {
                if (output_states[0] == CERTIFIED
                    && output_states[level] == SELECTED) {
                    input_states[0] = PENDING;
                }
                if (output_states[level] == CERTIFIED
                    && output_states[0] == SELECTED) {
                    input_states[level] = PENDING;
                }
            }
            return true;
        }

        if (output_states[0] != output_states[level]) {
            return false;
        }
        const std::size_t state = output_states[0];
        if (state != CERTIFIED) {
            return transition == 0;
        }
        if (transition == 0) {
            input_states[0] = CERTIFIED;
            input_states[level] = PENDING;
            return true;
        }
        if (transition == 1) {
            input_states[0] = PENDING;
            input_states[level] = CERTIFIED;
            return true;
        }
        if (transition == 2) {
            input_states[0] = CERTIFIED;
            input_states[level] = CERTIFIED;
            return true;
        }
        return false;
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
        if (merge_state == SELECTED) {
            first_state = SELECTED;
            second_state = SELECTED;
            return true;
        }
        if (merge_state == CERTIFIED) {
            first_state = CERTIFIED;
            second_state = PENDING;
            return true;
        }
        if (merge_state == PENDING) {
            first_state = PENDING;
            second_state = CERTIFIED;
            return true;
        }
        return false;
    }

    bool mergeInputStates(
            NetworKit::node left_vertex, NetworKit::node middle_vertex,
            NetworKit::node right_vertex,
            std::size_t output_left_state, std::size_t merge_state,
            std::size_t output_right_state,
            std::size_t &first_left_state,
            std::size_t &first_middle_state,
            std::size_t &second_middle_state,
            std::size_t &second_right_state) const {
        if (output_left_state >= stateCount()
            || merge_state >= stateCount()
            || output_right_state >= stateCount()) {
            return false;
        }
        const bool preserved_on_left = left_vertex == middle_vertex;
        const bool preserved_on_right = middle_vertex == right_vertex;
        if (!preserved_on_left && !preserved_on_right) {
            first_left_state = output_left_state;
            second_right_state = output_right_state;
            return splitMergeState(
                merge_state, first_middle_state, second_middle_state);
        }

        if (preserved_on_left && preserved_on_right) {
            return mergeFullyPreservedState(
                output_left_state, merge_state, output_right_state,
                first_left_state, first_middle_state,
                second_middle_state, second_right_state);
        }
        if (preserved_on_left) {
            second_right_state = output_right_state;
            return mergeLeftPreservedState(
                output_left_state, merge_state,
                first_left_state, first_middle_state,
                second_middle_state);
        }
        first_left_state = output_left_state;
        return mergeRightPreservedState(
            merge_state, output_right_state, first_middle_state,
            second_middle_state, second_right_state);
    }

    bool canonicalMergeOutputStates(
            NetworKit::node left_vertex, NetworKit::node right_vertex,
            std::size_t left_state, std::size_t right_state,
            std::size_t &output_left_state,
            std::size_t &output_right_state) const {
        if (left_state >= stateCount() || right_state >= stateCount()) {
            return false;
        }
        output_left_state = left_state;
        output_right_state = right_state;
        if (left_vertex != right_vertex) {
            return true;
        }

        const bool left_selected = left_state == SELECTED;
        const bool right_selected = right_state == SELECTED;
        if (left_selected != right_selected) {
            return false;
        }
        if (!left_selected
            && (left_state == CERTIFIED || right_state == CERTIFIED)) {
            output_left_state = CERTIFIED;
            output_right_state = CERTIFIED;
        }
        return true;
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
        return state == SELECTED || state == CERTIFIED;
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
            || output_states[0] != output_states[output_level]
            || !validStates(output_states)) {
            return false;
        }

        const std::size_t introduced_state = output_states[0];
        if (introduced_state == CERTIFIED
            && !introducedHasSelectedNeighbor(
                graph, vertex, input_boundary, output_states)) {
            return false;
        }

        input_states.resize(2 * input_level);
        for (std::size_t i = 0; i < input_level; ++i) {
            input_states[i] = precedingState(
                graph, vertex, input_boundary.left[i], introduced_state,
                output_states[i + 1]);
            input_states[input_level + i] = precedingState(
                graph, vertex, input_boundary.right[i], introduced_state,
                output_states[output_level + i + 1]);
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
                states.begin(), states.end(), [](std::size_t state) {
                    return state == SELECTED || state == CERTIFIED;
                });
    }

    bool better(Value first, Value second) const {
        return first < second;
    }

    void appendToSolution(
            NetworKit::node vertex, std::size_t state,
            Solution &solution) const {
        if (state == SELECTED) {
            solution.insert(vertex);
        }
    }

    void finalizeSolution(Solution &) const { }

 private:
    bool mergeLeftPreservedState(
            std::size_t output_state, std::size_t merge_state,
            std::size_t &first_left_state,
            std::size_t &first_middle_state,
            std::size_t &second_middle_state) const {
        if (output_state == SELECTED) {
            if (merge_state != SELECTED) {
                return false;
            }
            first_left_state = SELECTED;
            first_middle_state = SELECTED;
            second_middle_state = SELECTED;
            return true;
        }
        if (output_state == PENDING) {
            if (merge_state == SELECTED) {
                first_left_state = CERTIFIED;
                first_middle_state = CERTIFIED;
                second_middle_state = PENDING;
                return true;
            }
            if (merge_state == PENDING) {
                first_left_state = PENDING;
                first_middle_state = PENDING;
                second_middle_state = CERTIFIED;
                return true;
            }
            first_left_state = PENDING;
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            return merge_state == CERTIFIED;
        }
        if (merge_state == SELECTED) {
            first_left_state = CERTIFIED;
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            return true;
        }
        if (merge_state == CERTIFIED) {
            first_left_state = CERTIFIED;
            first_middle_state = CERTIFIED;
            second_middle_state = PENDING;
            return true;
        }
        first_left_state = PENDING;
        first_middle_state = PENDING;
        second_middle_state = CERTIFIED;
        return merge_state == PENDING;
    }

    bool mergeRightPreservedState(
            std::size_t merge_state, std::size_t output_state,
            std::size_t &first_middle_state,
            std::size_t &second_middle_state,
            std::size_t &second_right_state) const {
        if (output_state == SELECTED) {
            if (merge_state != SELECTED) {
                return false;
            }
            first_middle_state = SELECTED;
            second_middle_state = SELECTED;
            second_right_state = SELECTED;
            return true;
        }
        if (output_state == PENDING) {
            if (merge_state == SELECTED) {
                first_middle_state = PENDING;
                second_middle_state = CERTIFIED;
                second_right_state = CERTIFIED;
                return true;
            }
            if (merge_state == PENDING) {
                first_middle_state = CERTIFIED;
                second_middle_state = PENDING;
                second_right_state = PENDING;
                return true;
            }
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            second_right_state = PENDING;
            return merge_state == CERTIFIED;
        }
        if (merge_state == SELECTED) {
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            second_right_state = CERTIFIED;
            return true;
        }
        if (merge_state == CERTIFIED) {
            first_middle_state = CERTIFIED;
            second_middle_state = PENDING;
            second_right_state = PENDING;
            return true;
        }
        first_middle_state = PENDING;
        second_middle_state = CERTIFIED;
        second_right_state = CERTIFIED;
        return merge_state == PENDING;
    }

    bool mergeFullyPreservedState(
            std::size_t output_left_state, std::size_t merge_state,
            std::size_t output_right_state,
            std::size_t &first_left_state,
            std::size_t &first_middle_state,
            std::size_t &second_middle_state,
            std::size_t &second_right_state) const {
        const bool left_selected = output_left_state == SELECTED;
        const bool right_selected = output_right_state == SELECTED;
        if (left_selected != right_selected) {
            return false;
        }
        if (left_selected) {
            if (merge_state != SELECTED) {
                return false;
            }
            first_left_state = SELECTED;
            first_middle_state = SELECTED;
            second_middle_state = SELECTED;
            second_right_state = SELECTED;
            return true;
        }
        const bool requires_certificate =
            output_left_state == CERTIFIED
            || output_right_state == CERTIFIED;
        if (!requires_certificate) {
            if (merge_state == SELECTED) {
                first_left_state = CERTIFIED;
                first_middle_state = CERTIFIED;
                second_middle_state = PENDING;
                second_right_state = PENDING;
                return true;
            }
            if (merge_state == PENDING) {
                first_left_state = PENDING;
                first_middle_state = PENDING;
                second_middle_state = CERTIFIED;
                second_right_state = CERTIFIED;
                return true;
            }
            first_left_state = PENDING;
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            second_right_state = PENDING;
            return merge_state == CERTIFIED;
        }
        if (merge_state == SELECTED) {
            first_left_state = CERTIFIED;
            first_middle_state = PENDING;
            second_middle_state = PENDING;
            second_right_state = PENDING;
            return true;
        }
        if (merge_state == CERTIFIED) {
            first_left_state = CERTIFIED;
            first_middle_state = CERTIFIED;
            second_middle_state = PENDING;
            second_right_state = PENDING;
            return true;
        }
        if (merge_state == PENDING) {
            first_left_state = PENDING;
            first_middle_state = PENDING;
            second_middle_state = CERTIFIED;
            second_right_state = CERTIFIED;
            return true;
        }
        return false;
    }

    bool validStates(const std::vector<std::size_t> &states) const {
        return std::all_of(
            states.begin(), states.end(), [this](std::size_t state) {
                return state < stateCount();
            });
    }

    bool hasSelectedNeighbor(
            const NetworKit::Graph &graph, NetworKit::node vertex,
            const std::vector<NetworKit::node> &vertices,
            const std::vector<std::size_t> &states) const {
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (states[i] == SELECTED
                && graph.hasEdge(vertex, vertices[i])) {
                return true;
            }
        }
        return false;
    }

    bool introducedHasSelectedNeighbor(
            const NetworKit::Graph &graph, NetworKit::node vertex,
            const BakerSliceBoundary &input_boundary,
            const std::vector<std::size_t> &output_states) const {
        const std::size_t input_level = input_boundary.left.size();
        const std::size_t output_level = input_level + 1;
        for (std::size_t i = 0; i < input_level; ++i) {
            if (output_states[i + 1] == SELECTED
                && graph.hasEdge(vertex, input_boundary.left[i])) {
                return true;
            }
            if (output_states[output_level + i + 1] == SELECTED
                && graph.hasEdge(vertex, input_boundary.right[i])) {
                return true;
            }
        }
        return false;
    }

    std::size_t precedingState(
            const NetworKit::Graph &graph, NetworKit::node introduced,
            NetworKit::node boundary_vertex,
            std::size_t introduced_state,
            std::size_t output_state) const {
        if (output_state == CERTIFIED
            && introduced_state == SELECTED
            && graph.hasEdge(introduced, boundary_vertex)) {
            return PENDING;
        }
        return output_state;
    }
};

}  // namespace Koala
