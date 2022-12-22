/*
 * OddHoles.cpp
 *
 *  Created on: 05.11.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include <graph/GraphTools.hpp>
#include <traversal/BFS.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

#include "../other/oddHoles.h"

namespace Koala {

bool is_path(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &path) {
    for (unsigned i = 0; i < path.size() - 1; i++) {
        if (!graph.hasEdge(path[i], path[i + 1])) {
            return false;
        }
        for (unsigned j = i + 2; j < path.size(); j++) {
            if (graph.hasEdge(path[i], path[j])) {
                return false;
            }
        }
    }
    return true;
}

bool is_complete(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X, NetworKit::node v) {
    return std::none_of(X.begin(), X.end(), [&](auto i) {
        return v == i || !graph.hasEdge(v, i);
    });
}

std::vector<NetworKit::node> get_all_complete_vertices(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X) {
    std::vector<NetworKit::node> out;
    for (const auto &v : graph.nodeRange()) {
        if (is_complete(graph, X, v)) {
            out.push_back(v);
        }
    }
    return out;
}

bool PerfectGraphRecognition::containsOddHole(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(
        graph, std::numeric_limits<NetworKit::count>::max(), path,
        Koala::Traversal::PathInplaceMode::INDUCED_ODD_HOLE);
}

bool PerfectGraphRecognition::containsHole(const NetworKit::Graph &graph, NetworKit::count length) {
    if (length <= 3) {
        return false;
    }
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(
        graph, length, path, Koala::Traversal::PathInplaceMode::INDUCED_CYCLE);
}

bool PerfectGraphRecognition::containsT1(const NetworKit::Graph &graph) {
    //TODO(kturowski): temporary check
    Graph G(graph);
    assert(containsHole(graph, 5) == ::containsT1(G));
    if (containsHole(graph, 5)) {
        assert(containsOddHole(graph));
    }
    return containsHole(graph, 5);
}

bool PerfectGraphRecognition::containsT2(const NetworKit::Graph &graph) {
    //TODO(kturowski): temporary check
    Graph G(graph);
    for (const auto &v1 : graph.nodeRange()) {
        for (const auto &v2 : graph.neighborRange(v1)) {
            for (const auto &v3 : graph.neighborRange(v2)) {
                if (v3 == v1) {
                    continue;
                }
                for (const auto &v4 : graph.neighborRange(v3)) {
                    if (v4 == v1 || v4 == v2) {
                        continue;
                    }
                    if (!is_path(graph, std::vector<NetworKit::node>{v1, v2, v3, v4})) {
                        continue;
                    }
                    auto Y = get_all_complete_vertices(graph, {v1, v2, v4});
                    auto auxiliary_graph = NetworKit::GraphTools::subgraphFromNodes(
                        Koala::GraphTools::toComplement(graph), Y.begin(), Y.end());
                    NetworKit::ConnectedComponents auxiliary_components(auxiliary_graph);
                    auxiliary_components.run();
                    for (auto X : auxiliary_components.getComponents()) {
                        if (X.empty()) {
                            continue;
                        }
                        bool path = Koala::Traversal::BFS(graph, v1, v4, [&](auto v) {
                            if (v == v2 || v == v3) {
                                return false;
                            }
                            if (graph.hasEdge(v, v2) || graph.hasEdge(v, v3)) {
                                return false;
                            }
                            if (is_complete(graph, X, v)) {
                                return false;
                            }
                            return true;
                        });
                        if (path) {
                            assert(::containsT2(G));
                            return true;
                        }
                    }
                }
            }
        }
    }
    assert(!::containsT2(G));
    return false;
}

}
