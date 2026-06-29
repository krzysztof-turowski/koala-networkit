/*
 * ChordalGraphRecognition.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include "recognition/ChordalGraphRecognition.hpp"

#include <vector>

namespace Koala {

PerfectEliminationOrdering::PerfectEliminationOrdering(NetworKit::count n, NetworKit::count bound)
    : alpha(bound, 0), alpha_inv(n + 1, NetworKit::none) {
}

void PerfectEliminationOrdering::set(NetworKit::node v, NetworKit::index pos) {
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

const PerfectEliminationOrdering& ChordalGraphRecognition::getPEO() const {
    assureFinished();
    if (is_chordal != State::CHORDAL) {
        throw std::logic_error(
            "Cannot retrieve a Perfect Elimination Ordering: The graph is not chordal.");
    }
    return peo;
}

} /* namespace Koala */
