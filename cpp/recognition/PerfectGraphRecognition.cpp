/*
 * PerfectGraphRecognition.cpp
 *
 *  Created on: 11.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <list>
#include <set>

#include <iostream>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include <graph/GraphTools.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

#include "other/perfect.h"

namespace Koala {

PerfectGraphRecognition::PerfectGraphRecognition(
        const NetworKit::Graph &graph) : graph(std::make_optional(graph)), is_perfect(false) { }

bool PerfectGraphRecognition::isPerfect() const {
    assureFinished();
    return is_perfect;
}

void PerfectGraphRecognition::run() {
    is_perfect = isPerfectGraph(*graph);
    hasRun = true;
}

}  // namespace Koala
