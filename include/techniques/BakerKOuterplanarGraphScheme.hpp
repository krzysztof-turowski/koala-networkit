#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

enum class BakerFakeEdgeType {
    BRIDGE_COPY,
    CONNECTIVITY,
    TRIANGULATION
};

struct BakerFakeEdge {
    NetworKit::node u;
    NetworKit::node v;
    BakerFakeEdgeType type;
};

enum class BakerFaceTreeVertexType {
    FACE,
    EXTERIOR_EDGE,
    SINGLETON
};

struct BakerSliceBoundary {
    std::vector<NetworKit::node> left;
    std::vector<NetworKit::node> right;

    std::size_t size() const {
        return left.size() + right.size();
    }
};

/**
 * A vertex of one of Baker's ordered trees for a level component.
 */
struct BakerFaceTreeVertex {
    BakerFaceTreeVertexType type = BakerFaceTreeVertexType::FACE;
    NetworKit::node labelLeft = NetworKit::none;
    NetworKit::node labelRight = NetworKit::none;
    std::optional<std::size_t> parent;
    std::vector<std::size_t> children;
    std::optional<std::size_t> enclosedTree;
    std::size_t leftBoundaryIndex = 0;
    std::size_t rightBoundaryIndex = 0;
    std::size_t middleBoundaryIndex = 0;
    BakerSliceBoundary boundary;
    bool originalLabelEdge = false;
};

/**
 * The ordered, rooted face tree of one connected component at one Baker level.
 */
struct BakerFaceTree {
    NetworKit::count level = 0;
    std::vector<NetworKit::node> component;
    std::vector<BakerFaceTreeVertex> vertices;
    std::vector<std::size_t> leaves;
    std::size_t root = 0;
    std::optional<std::size_t> enclosingTree;
    std::optional<std::size_t> enclosingVertex;
};

struct BakerOperationStatistics {
    std::size_t recursiveCalls = 0;
    std::size_t adjustCalls = 0;
    std::size_t mergeCalls = 0;
    std::size_t contractCalls = 0;
    std::size_t createCalls = 0;
    std::size_t extendCalls = 0;
    std::size_t tableEntries = 0;
    std::size_t mergeTransitions = 0;
    std::size_t maximumTableEntries = 0;
    std::size_t maximumMergeTransitions = 0;
};

/**
 * Rotation system of a planar embedding, indexed by graph node identifier.
 */
struct BakerPlaneEmbedding {
    std::vector<std::vector<NetworKit::node>> rotation;
};

/**
 * Baker's forest of ordered level-component trees and their explicit slices.
 */
class BakerForest {
 public:
    std::size_t width() const {
        std::size_t result = 0;
        for (const auto &tree : trees) {
            for (const auto &vertex : tree.vertices) {
                result = std::max(result, vertex.boundary.size());
            }
        }
        return result;
    }

    NetworKit::count outerplanarity() const {
        return outerplanarity_;
    }

    bool hasTwoKBoundaryBound() const {
        for (const auto &tree : trees) {
            if (tree.level > outerplanarity_) {
                return false;
            }
            for (const auto &vertex : tree.vertices) {
                if (vertex.boundary.left.size() != tree.level
                    || vertex.boundary.right.size() != tree.level
                    || vertex.boundary.size() > 2 * outerplanarity_) {
                    return false;
                }
                for (std::size_t i = 0; i < tree.level; ++i) {
                    const auto expected_level = tree.level - i;
                    if (vertex.boundary.left[i] >= levels.size()
                        || vertex.boundary.right[i] >= levels.size()
                        || levels[vertex.boundary.left[i]] != expected_level
                        || levels[vertex.boundary.right[i]]
                            != expected_level) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    std::vector<BakerFaceTree> trees;
    std::vector<std::size_t> roots;
    std::vector<NetworKit::count> levels;
    BakerPlaneEmbedding embedding;
    std::vector<BakerFakeEdge> fakeEdges;
    BakerOperationStatistics statistics;

 private:
    NetworKit::count outerplanarity_ = 0;

    friend BakerForest buildBakerForest(const NetworKit::Graph &graph);
    friend BakerForest buildBakerForest(
        const NetworKit::Graph &graph,
        const std::vector<NetworKit::node> &outer_face_vertices);
    friend BakerForest buildBakerForest(
        const NetworKit::Graph &graph,
        const BakerPlaneEmbedding &embedding,
        const std::vector<NetworKit::node> &outer_face_vertices);
};

/**
 * Construct the ordered face trees, LB/RB maps, and slices from Baker's paper.
 */
BakerForest buildBakerForest(const NetworKit::Graph &graph);

/**
 * Construct a Baker forest using one prescribed outer face per component.
 *
 * The intersection of outer_face_vertices with each connected component must
 * be exactly the vertex set of a face in the computed planar embedding.
 */
BakerForest buildBakerForest(
    const NetworKit::Graph &graph,
    const std::vector<NetworKit::node> &outer_face_vertices);

/**
 * Construct a Baker forest from a prescribed embedding and outer faces.
 */
BakerForest buildBakerForest(
    const NetworKit::Graph &graph,
    const BakerPlaneEmbedding &embedding,
    const std::vector<NetworKit::node> &outer_face_vertices);

/**
 * General dynamic program over Baker's explicit k-outerplanar slices.
 *
 * Problem supplies the finite state semantics for each named table operation.
 * The scheme itself only constructs Baker's slices and enumerates their bounded
 * boundary assignments. Fake embedding edges are never treated as original
 * graph edges by a problem adapter.
 *
 * @tparam Problem Finite-state, additive boundary-table problem adapter.
 */
template <typename Problem>
class BakerKOuterplanarGraphScheme {
 public:
    using Value = typename Problem::Value;
    using Solution = typename Problem::Solution;

    BakerKOuterplanarGraphScheme() = default;

    explicit BakerKOuterplanarGraphScheme(const NetworKit::Graph &graph)
        : input_graph_(graph) { }

    explicit BakerKOuterplanarGraphScheme(Problem problem)
        : problem_(std::move(problem)) { }

    BakerKOuterplanarGraphScheme(const NetworKit::Graph &graph, Problem problem)
        : input_graph_(graph), problem_(std::move(problem)) { }

    /**
     * Solve the graph supplied to the constructor.
     *
     * @return A solution reconstructed by Problem.
     */
    Solution solve() {
        if (!input_graph_.has_value()) {
            throw std::logic_error("No graph was supplied to the Baker scheme");
        }
        return solve(*input_graph_);
    }

    /**
     * Construct Baker's embedding forest and solve the supplied graph.
     *
     * @param graph Undirected planar graph.
     * @return A solution reconstructed by Problem.
     */
    Solution solve(const NetworKit::Graph &graph) {
        baker_forest_ = buildBakerForest(graph);
        return solvePreparedForest(graph);
    }

    /**
     * Construct and solve using prescribed outer-face vertices.
     *
     * @param graph Undirected planar graph.
     * @param outer_face_vertices Union of one outer-face vertex set for each
     *        connected component.
     * @return A solution reconstructed by Problem.
     */
    Solution solve(
            const NetworKit::Graph &graph,
            const std::vector<NetworKit::node> &outer_face_vertices) {
        baker_forest_ = buildBakerForest(graph, outer_face_vertices);
        return solvePreparedForest(graph);
    }

    /**
     * Construct and solve using a prescribed embedding and outer faces.
     */
    Solution solve(
            const NetworKit::Graph &graph,
            const BakerPlaneEmbedding &embedding,
            const std::vector<NetworKit::node> &outer_face_vertices) {
        baker_forest_ = buildBakerForest(
            graph, embedding, outer_face_vertices);
        return solvePreparedForest(graph);
    }

 private:
    Solution solvePreparedForest(const NetworKit::Graph &graph) {
        if (!baker_forest_.hasTwoKBoundaryBound()) {
            throw std::logic_error("A Baker slice exceeded its formal 2k boundary");
        }

        graph_ = &graph;
        state_count_ = problem_.stateCount();
        if (state_count_ == 0) {
            throw std::logic_error("A Baker problem must define at least one state");
        }
        if (problem_.adjustmentCount() == 0
            || problem_.adjustmentCount() > state_count_) {
            throw std::logic_error(
                "A Baker adapter exceeded its bounded adjust transitions");
        }

        const std::size_t tree_count = baker_forest_.trees.size();
        tables_.clear();
        call_state_.clear();
        tables_.resize(tree_count);
        call_state_.resize(tree_count);
        std::size_t number_of_tree_vertices = 0;
        for (std::size_t tree = 0; tree < tree_count; ++tree) {
            const std::size_t size = baker_forest_.trees[tree].vertices.size();
            tables_[tree].resize(size);
            call_state_[tree].assign(size, 0);
            number_of_tree_vertices += size;
        }
        baker_forest_.statistics = {};

        Solution solution;
        for (const auto root_tree : baker_forest_.roots) {
            const auto &tree = baker_forest_.trees[root_tree];
            const auto table = computeTable(root_tree, tree.root);
            std::size_t best_code = invalid_index_;
            Value best_value = problem_.infeasibleValue();
            std::vector<std::size_t> states(2 * table->level, 0);
            for (std::size_t code = 0; code < table->values.size(); ++code) {
                const Value &candidate = table->values[code];
                if (problem_.isInfeasible(candidate)) {
                    continue;
                }
                decode(code, states);
                if (!problem_.isFinalAssignment(
                        table->boundary, states)) {
                    continue;
                }
                if (best_code == invalid_index_
                    || problem_.better(candidate, best_value)) {
                    best_code = code;
                    best_value = candidate;
                }
            }
            if (best_code == invalid_index_) {
                throw std::logic_error("The Baker problem has no feasible solution");
            }
            reconstruct(table, best_code, solution);
        }

        if (baker_forest_.statistics.recursiveCalls != number_of_tree_vertices) {
            throw std::logic_error("The Baker recursion did not visit every tree vertex once");
        }
        verifyComplexityBound(number_of_tree_vertices);
        problem_.finalizeSolution(solution);
        graph_ = nullptr;
        return solution;
    }

 public:
    const BakerForest &getBakerForest() const {
        return baker_forest_;
    }

 protected:
    struct Decision {
        std::size_t first = invalid_index_;
        std::size_t second = invalid_index_;
    };

    enum class TableOperation {
        DIRECT,
        ADJUST,
        MERGE,
        CONTRACT,
        CREATE,
        EXTEND
    };

    struct Table {
        NetworKit::count level = 0;
        BakerSliceBoundary boundary;
        std::vector<Value> values;
        std::vector<Decision> decisions;
        TableOperation operation = TableOperation::DIRECT;
        std::shared_ptr<Table> first;
        std::shared_ptr<Table> second;
        std::optional<NetworKit::node> introduced;
    };

    using TablePointer = std::shared_ptr<Table>;

    /**
     * Incorporate a face-tree label into a table.
     */
    TablePointer adjust(
            const TablePointer &input, NetworKit::node x, NetworKit::node y) {
        ++baker_forest_.statistics.adjustCalls;
        auto output = makeEmptyTable(
            input->level, input->boundary, TableOperation::ADJUST);
        output->first = input;

        std::vector<std::size_t> output_states(2 * input->level, 0);
        std::vector<std::size_t> input_states(2 * input->level, 0);
        const std::size_t adjustment_count = problem_.adjustmentCount();
        for (std::size_t output_code = 0;
             output_code < output->values.size(); ++output_code) {
            decode(output_code, output_states);
            for (std::size_t transition = 0;
                 transition < adjustment_count; ++transition) {
                if (!problem_.adjustInputStates(
                        *graph_, input->boundary, x, y, output_states,
                        transition, input_states)) {
                    continue;
                }
                const std::size_t input_code = encode(input_states);
                const Value &input_value = input->values[input_code];
                if (problem_.isInfeasible(input_value)) {
                    continue;
                }
                Value candidate = problem_.adjustValue(
                    input_value, x, y, output_states);
                if (problem_.isInfeasible(output->values[output_code])
                    || problem_.better(
                        candidate, output->values[output_code])) {
                    output->values[output_code] = std::move(candidate);
                    output->decisions[output_code].first = input_code;
                }
            }
        }
        return output;
    }

    /**
     * Merge two same-level slices sharing one complete boundary.
     */
    TablePointer merge(const TablePointer &left, const TablePointer &right) {
        ++baker_forest_.statistics.mergeCalls;
        if (left->level != right->level
            || left->boundary.right != right->boundary.left) {
            throw std::logic_error("Baker merge requires one identical boundary");
        }

        BakerSliceBoundary boundary{
            left->boundary.left, right->boundary.right};
        auto output = makeEmptyTable(
            left->level, std::move(boundary), TableOperation::MERGE);
        output->first = left;
        output->second = right;

        const std::size_t side_assignments = assignmentCount(left->level);
        std::vector<std::size_t> left_states(left->level, 0);
        std::vector<std::size_t> merge_states(left->level, 0);
        std::vector<std::size_t> output_left_states(left->level, 0);
        std::vector<std::size_t> output_right_states(left->level, 0);
        std::vector<std::size_t> first_left_states(left->level, 0);
        std::vector<std::size_t> first_middle_states(left->level, 0);
        std::vector<std::size_t> second_middle_states(left->level, 0);
        std::vector<std::size_t> second_right_states(left->level, 0);
        std::vector<std::size_t> right_states(left->level, 0);

        for (std::size_t left_code = 0;
             left_code < side_assignments; ++left_code) {
            decodeSide(left_code, left_states);
            for (std::size_t right_code = 0;
                 right_code < side_assignments; ++right_code) {
                decodeSide(right_code, right_states);
                // Canonicalization may identify two surviving occurrences
                // without adding a state digit to the merge enumeration.
                bool canonical = true;
                for (std::size_t i = 0; i < left->level; ++i) {
                    if (!problem_.canonicalMergeOutputStates(
                            left->boundary.left[i],
                            right->boundary.right[i], left_states[i],
                            right_states[i], output_left_states[i],
                            output_right_states[i])) {
                        canonical = false;
                        break;
                    }
                }
                if (!canonical) {
                    continue;
                }
                const std::size_t output_code = combineSideCodes(
                    encode(output_left_states), encode(output_right_states),
                    left->level);
                for (std::size_t middle_code = 0;
                     middle_code < side_assignments; ++middle_code) {
                    ++baker_forest_.statistics.mergeTransitions;
                    decodeSide(middle_code, merge_states);
                    bool feasible = true;
                    for (std::size_t i = 0; i < left->level; ++i) {
                        if (!problem_.mergeInputStates(
                                left->boundary.left[i],
                                left->boundary.right[i],
                                right->boundary.right[i], left_states[i],
                                merge_states[i], right_states[i],
                                first_left_states[i], first_middle_states[i],
                                second_middle_states[i],
                                second_right_states[i])) {
                            feasible = false;
                            break;
                        }
                    }
                    if (!feasible) {
                        continue;
                    }
                    const std::size_t first_left_code =
                        encode(first_left_states);
                    const std::size_t first_middle_code =
                        encode(first_middle_states);
                    const std::size_t second_middle_code =
                        encode(second_middle_states);
                    const std::size_t second_right_code =
                        encode(second_right_states);
                    const std::size_t first_code =
                        combineSideCodes(
                            first_left_code, first_middle_code, left->level);
                    const std::size_t second_code =
                        combineSideCodes(
                            second_middle_code, second_right_code,
                            left->level);
                    const Value &first_value = left->values[first_code];
                    const Value &second_value = right->values[second_code];
                    if (problem_.isInfeasible(first_value)
                        || problem_.isInfeasible(second_value)) {
                        continue;
                    }

                    // A Problem may use the bounded middle code to route an
                    // obligation, so correct duplicates from the actual child state.
                    Value candidate = problem_.mergeValue(
                        first_value, second_value,
                        left->boundary.right, first_middle_states);
                    if (problem_.isInfeasible(output->values[output_code])
                        || problem_.better(
                            candidate, output->values[output_code])) {
                        output->values[output_code] = std::move(candidate);
                        output->decisions[output_code] =
                            {first_code, second_code};
                    }
                }
            }
        }
        baker_forest_.statistics.maximumMergeTransitions = std::max(
            baker_forest_.statistics.maximumMergeTransitions,
            checkedPower(state_count_, 3 * left->level));
        return output;
    }

    /**
     * Remove the highest-level root vertex from both boundaries.
     */
    TablePointer contract(
            const TablePointer &input, const BakerSliceBoundary &boundary) {
        ++baker_forest_.statistics.contractCalls;
        if (input->level == 0
            || boundary.left.size() + 1 != input->boundary.left.size()
            || boundary.right.size() + 1 != input->boundary.right.size()
            || !std::equal(
                boundary.left.begin(), boundary.left.end(),
                input->boundary.left.begin() + 1)
            || !std::equal(
                boundary.right.begin(), boundary.right.end(),
                input->boundary.right.begin() + 1)) {
            throw std::logic_error("Invalid Baker contract boundaries");
        }

        const auto output_level =
            static_cast<NetworKit::count>(input->level - 1);
        auto output = makeEmptyTable(
            output_level, boundary, TableOperation::CONTRACT);
        output->first = input;
        std::vector<std::size_t> output_states(2 * output_level, 0);
        std::vector<std::size_t> input_states(2 * input->level, 0);

        for (std::size_t output_code = 0;
             output_code < output->values.size(); ++output_code) {
            decode(output_code, output_states);
            for (std::size_t top_state = 0;
                 top_state < state_count_; ++top_state) {
                if (!problem_.canContract(
                        input->boundary.left[0], top_state)) {
                    continue;
                }
                input_states[0] = top_state;
                input_states[input->level] = top_state;
                for (std::size_t i = 0; i < output_level; ++i) {
                    input_states[i + 1] = output_states[i];
                    input_states[input->level + i + 1] =
                        output_states[output_level + i];
                }
                const std::size_t input_code = encode(input_states);
                const Value &candidate = input->values[input_code];
                if (problem_.isInfeasible(candidate)) {
                    continue;
                }
                if (problem_.isInfeasible(output->values[output_code])
                    || problem_.better(
                        candidate, output->values[output_code])) {
                    output->values[output_code] = candidate;
                    output->decisions[output_code].first = input_code;
                }
            }
        }
        return output;
    }

    /**
     * Directly construct the small middle slice of rule S4.
     */
    TablePointer create(
            std::size_t tree_index, std::size_t vertex_index,
            std::size_t boundary_index) {
        ++baker_forest_.statistics.createCalls;
        const auto &tree = baker_forest_.trees[tree_index];
        const auto &vertex = tree.vertices[vertex_index];
        if (vertex.leftBoundaryIndex > boundary_index
            || boundary_index > vertex.rightBoundaryIndex) {
            throw std::logic_error("Baker create point lies outside LB/RB");
        }
        BakerSliceBoundary boundary = vertex.boundary;
        if (tree.level > 1) {
            if (!tree.enclosingTree.has_value()
                || !tree.enclosingVertex.has_value()) {
                throw std::logic_error(
                    "A Baker create operation has no enclosing face");
            }
            const auto &parent_tree =
                baker_forest_.trees[*tree.enclosingTree];
            const auto &parent_face =
                parent_tree.vertices[*tree.enclosingVertex];
            std::vector<NetworKit::node> lower_boundary;
            if (boundary_index == parent_face.children.size()) {
                lower_boundary = parent_tree.vertices[
                    parent_face.children.back()].boundary.right;
            } else {
                lower_boundary = parent_tree.vertices[
                    parent_face.children[boundary_index]].boundary.left;
            }
            boundary.left = {vertex.labelLeft};
            boundary.right = {vertex.labelRight};
            boundary.left.insert(
                boundary.left.end(), lower_boundary.begin(),
                lower_boundary.end());
            boundary.right.insert(
                boundary.right.end(), lower_boundary.begin(),
                lower_boundary.end());
        }
        return makeDirectTable(
            tree.level, boundary, TableOperation::CREATE);
    }

    /**
     * Add one next-level vertex to both boundaries of a table.
     */
    TablePointer extend(
            NetworKit::node vertex, const TablePointer &input) {
        ++baker_forest_.statistics.extendCalls;
        BakerSliceBoundary boundary;
        boundary.left.reserve(input->boundary.left.size() + 1);
        boundary.right.reserve(input->boundary.right.size() + 1);
        boundary.left.push_back(vertex);
        boundary.right.push_back(vertex);
        boundary.left.insert(
            boundary.left.end(), input->boundary.left.begin(),
            input->boundary.left.end());
        boundary.right.insert(
            boundary.right.end(), input->boundary.right.begin(),
            input->boundary.right.end());

        auto output = makeEmptyTable(
            input->level + 1, std::move(boundary), TableOperation::EXTEND);
        output->first = input;
        output->introduced = vertex;
        std::vector<std::size_t> output_states(2 * output->level, 0);
        std::vector<std::size_t> input_states(2 * input->level, 0);

        for (std::size_t code = 0; code < output->values.size(); ++code) {
            decode(code, output_states);
            if (!problem_.extendInputStates(
                    *graph_, vertex, input->boundary,
                    output_states, input_states)) {
                continue;
            }

            const std::size_t input_code = encode(input_states);
            const Value &input_value = input->values[input_code];
            if (problem_.isInfeasible(input_value)) {
                continue;
            }
            output->values[code] = problem_.extendValue(
                input_value, vertex, output_states[0]);
            output->decisions[code].first = input_code;
        }
        return output;
    }

    Problem problem_;

 private:
    static constexpr std::size_t invalid_index_ =
        std::numeric_limits<std::size_t>::max();

    std::size_t checkedPower(std::size_t base, std::size_t exponent) const {
        std::size_t result = 1;
        for (std::size_t i = 0; i < exponent; ++i) {
            if (result > std::numeric_limits<std::size_t>::max() / base) {
                throw std::overflow_error("A Baker table bound overflowed");
            }
            result *= base;
        }
        return result;
    }

    std::size_t checkedProduct(std::size_t first, std::size_t second) const {
        if (first != 0
            && second > std::numeric_limits<std::size_t>::max() / first) {
            throw std::overflow_error("A Baker operation bound overflowed");
        }
        return first * second;
    }

    std::size_t assignmentCount(std::size_t number_of_vertices) const {
        return checkedPower(state_count_, number_of_vertices);
    }

    void decode(std::size_t code, std::vector<std::size_t> &states) const {
        for (auto &state : states) {
            state = code % state_count_;
            code /= state_count_;
        }
    }

    void decodeSide(
            std::size_t code, std::vector<std::size_t> &states) const {
        decode(code, states);
    }

    std::size_t encode(const std::vector<std::size_t> &states) const {
        std::size_t result = 0;
        std::size_t multiplier = 1;
        for (const auto state : states) {
            if (state >= state_count_) {
                throw std::logic_error(
                    "A Baker adapter produced an unknown boundary state");
            }
            result += state * multiplier;
            multiplier *= state_count_;
        }
        return result;
    }

    std::size_t combineSideCodes(
            std::size_t left_code, std::size_t right_code,
            std::size_t level) const {
        return left_code + right_code * assignmentCount(level);
    }

    TablePointer makeEmptyTable(
            NetworKit::count level, BakerSliceBoundary boundary,
            TableOperation operation) {
        if (boundary.left.size() != level
            || boundary.right.size() != level
            || 2 * level > 2 * baker_forest_.outerplanarity()) {
            throw std::logic_error("A Baker table violated the 2k boundary");
        }
        const std::size_t table_size = assignmentCount(2 * level);
        const std::size_t maximum_size = assignmentCount(
            2 * baker_forest_.outerplanarity());
        if (table_size > maximum_size) {
            throw std::logic_error("A Baker table exceeded r^(2k) entries");
        }

        auto table = std::make_shared<Table>();
        table->level = level;
        table->boundary = std::move(boundary);
        table->operation = operation;
        table->values.assign(table_size, problem_.infeasibleValue());
        table->decisions.resize(table_size);
        baker_forest_.statistics.tableEntries += table_size;
        baker_forest_.statistics.maximumTableEntries = std::max(
            baker_forest_.statistics.maximumTableEntries, table_size);
        return table;
    }

    TablePointer makeDirectTable(
            NetworKit::count level, const BakerSliceBoundary &boundary,
            TableOperation operation) {
        auto table = makeEmptyTable(level, boundary, operation);
        std::vector<std::size_t> states(2 * level, 0);

        for (std::size_t code = 0; code < table->values.size(); ++code) {
            decode(code, states);
            const auto value = problem_.directValue(
                *graph_, boundary, states);
            if (value.has_value()) {
                table->values[code] = *value;
            }
        }
        return table;
    }

    TablePointer computeTable(
            std::size_t tree_index, std::size_t vertex_index) {
        if (call_state_[tree_index][vertex_index] != 0) {
            throw std::logic_error(
                "A Baker face-tree vertex was requested more than once");
        }
        call_state_[tree_index][vertex_index] = 1;
        ++baker_forest_.statistics.recursiveCalls;

        const auto &tree = baker_forest_.trees[tree_index];
        const auto &vertex = tree.vertices[vertex_index];
        TablePointer result;

        if (vertex.type == BakerFaceTreeVertexType::FACE) {
            if (vertex.enclosedTree.has_value()) {
                const auto child_tree = *vertex.enclosedTree;
                const auto child_root = baker_forest_.trees[child_tree].root;
                result = contract(
                    computeTable(child_tree, child_root), vertex.boundary);
                result = adjust(
                    result, vertex.labelLeft, vertex.labelRight);
            } else {
                if (vertex.children.empty()) {
                    throw std::logic_error(
                        "A Baker face-tree face has no children");
                }
                result = computeTable(tree_index, vertex.children.front());
                for (std::size_t i = 1; i < vertex.children.size(); ++i) {
                    result = merge(
                        result, computeTable(tree_index, vertex.children[i]));
                }
                result = adjust(
                    result, vertex.labelLeft, vertex.labelRight);
            }
        } else if (tree.level == 1) {
            result = makeDirectTable(
                tree.level, vertex.boundary, TableOperation::DIRECT);
        } else {
            if (!tree.enclosingTree.has_value()
                || !tree.enclosingVertex.has_value()) {
                throw std::logic_error(
                    "A non-exterior Baker tree has no enclosing face");
            }
            const auto parent_tree = *tree.enclosingTree;
            const auto parent_vertex = *tree.enclosingVertex;
            const auto &face =
                baker_forest_.trees[parent_tree].vertices[parent_vertex];
            const std::size_t left = vertex.leftBoundaryIndex;
            const std::size_t right = vertex.rightBoundaryIndex;
            if (left == right) {
                result = create(tree_index, vertex_index, left);
            } else {
                const std::size_t middle =
                    vertex.middleBoundaryIndex;
                if (middle < left || middle > right) {
                    throw std::logic_error(
                        "A Baker middle point lies outside LB/RB");
                }
                result = create(tree_index, vertex_index, middle);
                for (std::size_t index = middle; index > left; --index) {
                    const auto child = face.children[index - 1];
                    result = merge(
                        extend(
                            vertex.labelLeft,
                            computeTable(parent_tree, child)),
                        result);
                }
                for (std::size_t index = middle;
                     index < right; ++index) {
                    const auto child = face.children[index];
                    result = merge(
                        result,
                        extend(
                            vertex.labelRight,
                            computeTable(parent_tree, child)));
                }
            }
        }

        if (result->boundary.left != vertex.boundary.left
            || result->boundary.right != vertex.boundary.right) {
            throw std::logic_error(
                "A Baker table does not match its materialized slice");
        }
        tables_[tree_index][vertex_index] = result;
        call_state_[tree_index][vertex_index] = 2;
        return result;
    }

    void reconstruct(
            const TablePointer &table, std::size_t code,
            Solution &solution) const {
        if (code >= table->values.size()
            || problem_.isInfeasible(table->values[code])) {
            throw std::logic_error("Cannot reconstruct an infeasible Baker entry");
        }
        const Decision &decision = table->decisions[code];
        switch (table->operation) {
            case TableOperation::DIRECT:
            case TableOperation::CREATE: {
                std::vector<std::size_t> states(2 * table->level, 0);
                decode(code, states);
                for (std::size_t i = 0; i < table->level; ++i) {
                    problem_.appendToSolution(
                        table->boundary.left[i], states[i], solution);
                    problem_.appendToSolution(
                        table->boundary.right[i],
                        states[table->level + i], solution);
                }
                break;
            }
            case TableOperation::ADJUST:
            case TableOperation::CONTRACT:
                reconstruct(table->first, decision.first, solution);
                break;
            case TableOperation::MERGE:
                reconstruct(table->first, decision.first, solution);
                reconstruct(table->second, decision.second, solution);
                break;
            case TableOperation::EXTEND: {
                reconstruct(table->first, decision.first, solution);
                std::vector<std::size_t> states(2 * table->level, 0);
                decode(code, states);
                problem_.appendToSolution(
                    *table->introduced, states[0], solution);
                break;
            }
            default:
                throw std::logic_error("Unknown Baker table operation");
        }
    }

    void verifyComplexityBound(std::size_t number_of_tree_vertices) {
        const std::size_t k = baker_forest_.outerplanarity();
        const std::size_t table_bound = checkedPower(state_count_, 2 * k);
        const std::size_t merge_bound = checkedPower(state_count_, 3 * k);
        const auto &statistics = baker_forest_.statistics;
        const std::size_t operation_calls =
            statistics.adjustCalls + statistics.mergeCalls
            + statistics.contractCalls + statistics.createCalls
            + statistics.extendCalls;
        const std::size_t table_allocation_bound = checkedProduct(
            number_of_tree_vertices + operation_calls, table_bound);
        const std::size_t merge_transition_bound = checkedProduct(
            statistics.mergeCalls, merge_bound);
        if (statistics.maximumTableEntries > table_bound
            || statistics.maximumMergeTransitions > merge_bound
            || statistics.tableEntries > table_allocation_bound
            || statistics.mergeTransitions > merge_transition_bound
            || statistics.mergeCalls > number_of_tree_vertices
            || statistics.adjustCalls > number_of_tree_vertices
            || statistics.contractCalls > number_of_tree_vertices
            || statistics.createCalls > number_of_tree_vertices
            || statistics.extendCalls > number_of_tree_vertices) {
            throw std::logic_error(
                "The Baker dynamic program exceeded its formal complexity bound");
        }
    }

    std::optional<NetworKit::Graph> input_graph_;
    BakerForest baker_forest_;
    const NetworKit::Graph *graph_ = nullptr;
    std::size_t state_count_ = 0;
    std::vector<std::vector<TablePointer>> tables_;
    std::vector<std::vector<unsigned char>> call_state_;
};

}  // namespace Koala
