#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

#include "techniques/BakerKOuterplanarGraphScheme.hpp"

namespace Koala {

/**
 * Statistics for one execution of Baker's planar approximation scheme.
 */
struct BakerPlanarApproximationStatistics {
    NetworKit::count parameter = 0;
    NetworKit::count promisedOuterplanarity = 0;
    NetworKit::count maximumSubproblemOuterplanarity = 0;
    std::size_t residueCandidates = 0;
    std::size_t exactSubproblems = 0;
};

/**
 * Baker's planar layer approximation scheme.
 *
 * Maximization problems delete one residue modulo k + 1 and solve the
 * remaining k-outerplanar components. Minimization problems solve overlapping
 * strips of levels jk + i through (j + 1)k + i, inclusive, and unite their
 * solutions. Problem supplies the finite-state exact-scheme interface and the
 * optimization direction; common solution comparison, union, and node-ID
 * translation are implemented here.
 *
 * @tparam ExactScheme Exact bounded-outerplanarity scheme template.
 * @tparam Problem Exact and approximation problem adapter.
 */
template <template <typename> typename ExactScheme, typename Problem>
class BakerPlanarApproximationScheme {
 public:
    using Solution = typename Problem::Solution;

    const BakerPlanarApproximationStatistics &
    getBakerApproximationStatistics() const {
        return statistics_;
    }

    double getEpsilon() const {
        return epsilon_;
    }

 protected:
    BakerPlanarApproximationScheme(
            const NetworKit::Graph &graph, double epsilon)
        : graph_(graph), epsilon_(epsilon) {
        if (!std::isfinite(epsilon_) || epsilon_ <= 0.0) {
            throw std::invalid_argument("epsilon must be finite and positive");
        }
    }

    Solution solve() {
        statistics_ = {};
        if (graph_.isDirected()) {
            throw std::invalid_argument(
                "BakerPlanarApproximationScheme requires an undirected graph");
        }
        if (graph_.numberOfNodes() == 0) {
            return initialSolution();
        }

        const BakerForest level_forest = buildBakerForest(graph_);
        const auto &levels = level_forest.levels;
        const auto &embedding = level_forest.embedding;
        const NetworKit::count parameter = approximationParameter();
        statistics_.parameter = parameter;
        if (problem_.isMaximization()) {
            statistics_.promisedOuterplanarity = parameter;
            return solveMaximization(levels, embedding, parameter);
        }
        statistics_.promisedOuterplanarity = parameter + 1;
        return solveMinimization(levels, embedding, parameter);
    }

 private:
    struct CompactSubgraph {
        NetworKit::Graph graph;
        std::vector<NetworKit::node> localToGlobal;
        BakerPlaneEmbedding embedding;
        std::vector<NetworKit::node> outerFaceVertices;
    };

    NetworKit::count approximationParameter() const {
        const NetworKit::count number_of_nodes = graph_.numberOfNodes();
        const double reciprocal = 1.0 / epsilon_;
        if (!std::isfinite(reciprocal)
            || reciprocal >= static_cast<double>(number_of_nodes)) {
            return number_of_nodes;
        }
        return std::max<NetworKit::count>(
            1, static_cast<NetworKit::count>(std::ceil(reciprocal)));
    }

    CompactSubgraph compactInducedSubgraph(
            const std::vector<NetworKit::node> &vertices,
            const std::vector<NetworKit::count> &levels,
            const BakerPlaneEmbedding &embedding) const {
        CompactSubgraph result{
            NetworKit::Graph(vertices.size(), graph_.isWeighted(), false),
            vertices,
            BakerPlaneEmbedding{},
            {}};
        std::unordered_map<NetworKit::node, NetworKit::node> global_to_local;
        global_to_local.reserve(vertices.size());
        for (NetworKit::node local = 0; local < vertices.size(); ++local) {
            global_to_local.emplace(vertices[local], local);
        }

        for (NetworKit::node local_u = 0;
             local_u < vertices.size(); ++local_u) {
            const NetworKit::node global_u = vertices[local_u];
            graph_.forNeighborsOf(
                global_u,
                [&](NetworKit::node global_v, NetworKit::edgeweight weight) {
                    const auto local_v = global_to_local.find(global_v);
                    if (local_v == global_to_local.end()
                        || local_u > local_v->second) {
                        return;
                    }
                    result.graph.addEdge(local_u, local_v->second, weight);
                });
        }

        result.embedding.rotation.resize(vertices.size());
        for (NetworKit::node local_u = 0;
             local_u < vertices.size(); ++local_u) {
            const auto global_u = vertices[local_u];
            for (const auto global_v : embedding.rotation[global_u]) {
                const auto local_v = global_to_local.find(global_v);
                if (local_v != global_to_local.end()) {
                    result.embedding.rotation[local_u].push_back(
                        local_v->second);
                }
            }
        }

        std::vector<unsigned char> seen(vertices.size(), 0);
        for (NetworKit::node start = 0; start < vertices.size(); ++start) {
            if (seen[start]) {
                continue;
            }
            std::vector<NetworKit::node> component;
            NetworKit::count first_level = levels[vertices[start]];
            seen[start] = 1;
            component.push_back(start);
            for (std::size_t index = 0;
                 index < component.size(); ++index) {
                const auto u = component[index];
                first_level = std::min(first_level, levels[vertices[u]]);
                result.graph.forNeighborsOf(u, [&](NetworKit::node v) {
                    if (!seen[v]) {
                        seen[v] = 1;
                        component.push_back(v);
                    }
                });
            }
            for (const auto u : component) {
                if (levels[vertices[u]] == first_level) {
                    result.outerFaceVertices.push_back(u);
                }
            }
        }
        return result;
    }

    Solution solveInducedSubgraph(
            const std::vector<NetworKit::node> &vertices,
            const std::vector<NetworKit::count> &levels,
            const BakerPlaneEmbedding &embedding,
            NetworKit::count promised_outerplanarity) {
        if (vertices.empty()) {
            return initialSolution();
        }
        CompactSubgraph subgraph = compactInducedSubgraph(
            vertices, levels, embedding);
        ExactScheme<Problem> exact_scheme;
        Solution local_solution = exact_scheme.solve(
            subgraph.graph, subgraph.embedding,
            subgraph.outerFaceVertices);
        const NetworKit::count actual_outerplanarity =
            exact_scheme.getBakerForest().outerplanarity();
        if (actual_outerplanarity > promised_outerplanarity) {
            throw std::logic_error(
                "An exact Baker subproblem exceeded its promised layer width");
        }
        statistics_.maximumSubproblemOuterplanarity = std::max(
            statistics_.maximumSubproblemOuterplanarity,
            actual_outerplanarity);
        ++statistics_.exactSubproblems;
        return translateSolution(local_solution, subgraph.localToGlobal);
    }

    Solution solveMaximization(
            const std::vector<NetworKit::count> &levels,
            const BakerPlaneEmbedding &embedding,
            NetworKit::count parameter) {
        const NetworKit::count period = parameter + 1;
        std::optional<Solution> best;
        for (NetworKit::count residue = 0; residue < period; ++residue) {
            std::vector<NetworKit::node> vertices;
            vertices.reserve(graph_.numberOfNodes());
            graph_.forNodes([&](NetworKit::node u) {
                if (levels[u] % period != residue) {
                    vertices.push_back(u);
                }
            });
            Solution candidate = solveInducedSubgraph(
                vertices, levels, embedding, parameter);
            if (!best.has_value()
                || betterSolution(candidate, *best)) {
                best = std::move(candidate);
            }
            ++statistics_.residueCandidates;
        }
        return std::move(*best);
    }

    std::vector<std::vector<NetworKit::node>> levelBuckets(
            const std::vector<NetworKit::count> &levels,
            NetworKit::count &last_level) const {
        last_level = 0;
        graph_.forNodes([&](NetworKit::node u) {
            last_level = std::max(last_level, levels[u]);
        });
        std::vector<std::vector<NetworKit::node>> buckets(last_level + 1);
        graph_.forNodes([&](NetworKit::node u) {
            buckets[levels[u]].push_back(u);
        });
        return buckets;
    }

    std::vector<NetworKit::node> verticesInLevels(
            const std::vector<std::vector<NetworKit::node>> &buckets,
            NetworKit::count first, NetworKit::count last) const {
        std::vector<NetworKit::node> vertices;
        for (NetworKit::count level = first; level <= last; ++level) {
            vertices.insert(
                vertices.end(), buckets[level].begin(), buckets[level].end());
        }
        return vertices;
    }

    Solution solveMinimization(
            const std::vector<NetworKit::count> &levels,
            const BakerPlaneEmbedding &embedding,
            NetworKit::count parameter) {
        NetworKit::count last_level = 0;
        const auto buckets = levelBuckets(levels, last_level);
        std::optional<Solution> best;
        for (NetworKit::count residue = 0;
             residue < parameter; ++residue) {
            NetworKit::count last_band_level =
                residue == 0 ? parameter : residue;
            if (last_band_level == 1) {
                last_band_level += parameter;
            }

            Solution candidate = initialSolution();
            for (NetworKit::count first_level = 1;
                 first_level <= last_level;) {
                const NetworKit::count capped_last =
                    std::min(last_band_level, last_level);
                const auto vertices = verticesInLevels(
                    buckets, first_level, capped_last);
                Solution band_solution = solveInducedSubgraph(
                    vertices, levels, embedding, parameter + 1);
                mergeSolutions(candidate, band_solution);
                if (last_band_level > last_level) {
                    break;
                }
                first_level = last_band_level;
                last_band_level += parameter;
            }
            if (!best.has_value()
                || betterSolution(candidate, *best)) {
                best = std::move(candidate);
            }
            ++statistics_.residueCandidates;
        }
        return std::move(*best);
    }

    Solution initialSolution() const {
        return {};
    }

    bool betterSolution(
            const Solution &first, const Solution &second) const {
        return problem_.isMaximization()
            ? first.size() > second.size()
            : first.size() < second.size();
    }

    void mergeSolutions(
            Solution &target, const Solution &source) const {
        std::copy(
            source.begin(), source.end(),
            std::inserter(target, target.end()));
        problem_.finalizeSolution(target);
    }

    Solution translateSolution(
            const Solution &solution,
            const std::vector<NetworKit::node> &local_to_global) const {
        Solution result;
        auto output = std::inserter(result, result.end());
        for (const auto local : solution) {
            if (local >= local_to_global.size()) {
                throw std::logic_error(
                    "An exact Baker subproblem returned an unknown node");
            }
            *output++ = local_to_global[local];
        }
        problem_.finalizeSolution(result);
        return result;
    }

    const NetworKit::Graph &graph_;
    double epsilon_;
    Problem problem_;
    BakerPlanarApproximationStatistics statistics_;
};

}  // namespace Koala
