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

#include "other/perfect.h"

namespace Koala {

PerfectGraphRecognition::PerfectGraphRecognition(const NetworKit::Graph &graph)
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
    is_perfect = containsSimpleProhibited(*graph);
    if (is_perfect != State::UNKNOWN) {
        hasRun = true;
        return;
    }
    auto graph_complement = Koala::GraphTools::toComplement(*graph);
    is_perfect = containsSimpleProhibited(graph_complement);
    if (is_perfect != State::UNKNOWN) {
        hasRun = true;
        return;
    }

    is_perfect = checkPerfectGraph(*graph);
    hasRun = true;
}

void PerfectGraphRecognition::check() const {
    assureFinished();
    auto graph_complement = Koala::GraphTools::toComplement(*graph);
    assert((!containsOddHole(*graph) && !containsOddHole(graph_complement)) == isPerfect());
    //TODO(kturowski): temporary check
    assert(isPerfectGraphNaive(*graph, false) == isPerfect());
}

static PerfectGraphRecognition::State PerfectGraphRecognition::containsSimpleProhibited(
        const NetworKit::Graph &graph) {
    if (containsJewel(graph)) {
      return State::HAS_JEWEL;
    }
    if (containsPyramid(graph)) {
      return State::HAS_PYRAMID;
    }
    if (containsT1(graph)) {
      return State::HAS_T1;
    }
    if (containsT2(graph)) {
      return State::HAS_T2;
    }
    return State::UNKNOWN;
}

}  // namespace Koala
