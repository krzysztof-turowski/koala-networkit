/*
 * CographRecognition.cpp
 *
 *  Created on: 29.10.2023
 *      Author: fixikmila
 */

#include <graph/GraphTools.hpp>

#include "recognition/CographRecognition.hpp"

namespace Koala {

CographRecognition::CographRecognition(const NetworKit::Graph &graph)
    : graph(graph), is_cograph(State::UNKNOWN) { }

bool CographRecognition::isCograph() const {
    assureFinished();
    return is_cograph == State::COGRAPH;
}

CographRecognition::State CographRecognition::getState() const {
    assureFinished();
    return is_cograph;
}

bool check_path(
        const NetworKit::Graph &G,
        NetworKit::node x, NetworKit::node y, NetworKit::node u, NetworKit::node v) {
    return G.hasEdge(y, u) && !G.hasEdge(x, u) && !G.hasEdge(x, v) && !G.hasEdge(y, v);
}

void CographRecognition::check() const {
    assureFinished();
    for (const auto &[x, y] : graph.edgeRange()) {
        for (const auto &[u, v] : graph.edgeRange()) {
            if (x == u || x == v || y == u || y == v) {
                continue;
            }
            assert(!check_path(graph, x, y, u, v) && !check_path(graph, x, y, v, u)
                && !check_path(graph, y, x, u, v) && !check_path(graph, y, x, v, u));
        }
    }
}

}  // namespace Koala
