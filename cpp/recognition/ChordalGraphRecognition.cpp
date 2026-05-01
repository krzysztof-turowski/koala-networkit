/*
 * ChordalGraphRecognition.cpp
 *
 *  Created on:
 *      Author:
 */

#include "recognition/ChordalGraphRecognition.hpp"

#include <algorithm>
#include <vector>

namespace Koala {

PerfectEliminationOrdering::PerfectEliminationOrdering(NetworKit::count n, NetworKit::count bound)
    : alpha(bound, 0), alpha_inv(n + 1, NetworKit::none) {
}

void PerfectEliminationOrdering::set(NetworKit::node v, NetworKit::count pos) {
    alpha[v] = pos;
    alpha_inv[pos] = v;
}

ChordalGraphRecognition::ChordalGraphRecognition(const NetworKit::Graph& graph)
    : graph(graph), is_chordal(State::UNKNOWN) {
}

bool ChordalGraphRecognition::isChordal() const {
    assureFinished();
    return is_chordal == State::CHORDAL;
}

ChordalGraphRecognition::State ChordalGraphRecognition::getState() const {
    assureFinished();
    return is_chordal;
}

void ChordalGraphRecognition::check() const {
    assureFinished();
    // TODO
}

} /* namespace Koala */
