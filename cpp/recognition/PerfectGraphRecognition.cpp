/*
 * PerfectGraphRecognition.cpp
 *
 *  Created on: 11.11.2021
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <list>
#include <set>

#include <iostream>

#include <graph/GraphTools.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

namespace Koala {

PerfectGraphRecognition::PerfectGraphRecognition(NetworKit::Graph &graph)
    : graph(std::make_optional(graph)), is_perfect(State::UNKNOWN) { }

bool PerfectGraphRecognition::isPerfect() const {
    assureFinished();
    return is_perfect == State::PERFECT;
}

PerfectGraphRecognition::State PerfectGraphRecognition::getState() const {
    assureFinished();
    return is_perfect;
}

void PerfectGraphRecognition::run() {
    hasRun = true;
    if (graph->numberOfNodes() <= 4) {
        is_perfect = State::PERFECT;
        return;
    }
    is_perfect = contains_simple_prohibited(*graph);
    if (is_perfect != State::UNKNOWN) {
        return;
    }
    auto graph_complement = Koala::GraphTools::toComplement(*graph);
    is_perfect = contains_simple_prohibited(graph_complement);
    if (is_perfect != State::UNKNOWN) {
        return;
    }
    if (contains_near_cleaner_odd_hole(*graph)) {
        is_perfect = State::HAS_NEAR_CLEANER_ODD_HOLE;
        return;
    }
    if (contains_near_cleaner_odd_hole(graph_complement)) {
        is_perfect = State::HAS_NEAR_CLEANER_ODD_HOLE;
        return;
    }
    is_perfect = State::PERFECT;
}

void PerfectGraphRecognition::check() const {
    assureFinished();
    auto graph_complement = Koala::GraphTools::toComplement(*graph);
    assert((!contains_odd_hole(*graph) && !contains_odd_hole(graph_complement)) == isPerfect());
}

bool PerfectGraphRecognition::isComplete(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X, NetworKit::node v) {
    return std::none_of(X.begin(), X.end(), [&](auto i) {
        return v == i || !graph.hasEdge(v, i);
    });
}

std::vector<NetworKit::node> PerfectGraphRecognition::getAllCompleteVertices(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X) {
    std::vector<NetworKit::node> out;
    for (const auto &v : graph.nodeRange()) {
        if (isComplete(graph, X, v)) {
            out.push_back(v);
        }
    }
    return out;
}

std::vector<std::vector<NetworKit::node>> PerfectGraphRecognition::getAuxiliaryComponents(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &V) {
    auto Y = getAllCompleteVertices(graph, V);
    auto auxiliary_graph = NetworKit::GraphTools::subgraphFromNodes(
        Koala::GraphTools::toComplement(graph), Y.begin(), Y.end());
    NetworKit::ConnectedComponents auxiliary_components(auxiliary_graph);
    auxiliary_components.run();
    return auxiliary_components.getComponents();
}

PerfectGraphRecognition::State PerfectGraphRecognition::contains_simple_prohibited(
        const NetworKit::Graph &graph) {
    if (contains_jewel(graph)) {
        return State::HAS_JEWEL;
    }
    if (contains_pyramid(graph)) {
        return State::HAS_PYRAMID;
    }
    if (contains_t1(graph)) {
        return State::HAS_T1;
    }
    if (contains_t2(graph)) {
        return State::HAS_T2;
    }
    if (contains_t3(graph)) {
        return State::HAS_T3;
    }
    return State::UNKNOWN;
}

}  // namespace Koala
