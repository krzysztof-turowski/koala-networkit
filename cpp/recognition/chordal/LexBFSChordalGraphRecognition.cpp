/*
 * LexBFSChordalGraphRecognition.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include <algorithm>
#include <vector>

#include "recognition/ChordalGraphRecognition.hpp"

namespace Koala {

LexBFSChordalGraphRecognition::SetQueue::SetQueue(const NetworKit::Graph& graph) {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();

    sets.resize(n + 1);
    vertices.resize(bound);

    free_sets.reserve(n - 1);
    for (NetworKit::index i = n; i >= 2; i--) free_sets.push_back(i);

    first_set = INITIAL_SET;

    for (auto v : graph.nodeRange()) addVertex(v, INITIAL_SET);
}

void LexBFSChordalGraphRecognition::SetQueue::removeVertex(
    NetworKit::node v, NetworKit::index set) {
    NetworKit::node prev = vertices[v].prev;
    NetworKit::node next = vertices[v].next;

    if (prev != NetworKit::none) vertices[prev].next = next;
    else sets[set].head = next;

    if (next != NetworKit::none) vertices[next].prev = prev;

    vertices[v].prev = NetworKit::none;
    vertices[v].next = NetworKit::none;
}

void LexBFSChordalGraphRecognition::SetQueue::addVertex(NetworKit::node v, NetworKit::index set) {
    vertices[v].set_id = set;
    vertices[v].prev = NetworKit::none;

    NetworKit::node h = sets[set].head;
    vertices[v].next = h;

    if (h != NetworKit::none) vertices[h].prev = v;
    sets[set].head = v;
}

void LexBFSChordalGraphRecognition::SetQueue::removeSet(NetworKit::index set) {
    NetworKit::index prev = sets[set].prev;
    NetworKit::index next = sets[set].next;

    if (prev != NO_SET) sets[prev].next = next;
    else first_set = next;

    if (next != NO_SET) sets[next].prev = prev;

    free_sets.push_back(set);
}

NetworKit::node LexBFSChordalGraphRecognition::SetQueue::popVertex() {
    if (first_set == NO_SET) return NetworKit::none;

    NetworKit::node v = sets[first_set].head;

    removeVertex(v, first_set);
    vertices[v].set_id = PROCESSED_SET;

    if (sets[first_set].head == NetworKit::none) removeSet(first_set);

    return v;
}

void LexBFSChordalGraphRecognition::SetQueue::moveToNewSet(NetworKit::node w) {
    NetworKit::index s = vertices[w].set_id;
    if (s == PROCESSED_SET) return;

    NetworKit::index s_prime = sets[s].split_into;
    bool is_new_split = (s_prime == NO_SET);

    if (is_new_split) {
        s_prime = free_sets.back();
        free_sets.pop_back();

        NetworKit::index p = sets[s].prev;
        sets[s_prime].prev = p;
        sets[s_prime].next = s;
        sets[s].prev = s_prime;

        if (p != NO_SET) sets[p].next = s_prime;
        else first_set = s_prime;

        sets[s].split_into = s_prime;
        active_splits.push_back(s);
    }

    removeVertex(w, s);
    addVertex(w, s_prime);

    if (sets[s].head == NetworKit::none) removeSet(s);
}

void LexBFSChordalGraphRecognition::SetQueue::removeEmptySets() {
    for (NetworKit::index s : active_splits) {
        sets[s].split_into = NO_SET;
    }
    active_splits.clear();
}

bool LexBFSChordalGraphRecognition::SetQueue::isProcessed(NetworKit::node w) const {
    return vertices[w].set_id == PROCESSED_SET;
}

PerfectEliminationOrdering LexBFSChordalGraphRecognition::getLexBFSOrdering() const {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();

    PerfectEliminationOrdering ordering(n, bound);

    SetQueue queue(graph);

    for (NetworKit::index i = n; i >= 1; i--) {
        NetworKit::node v = queue.popVertex();

        ordering.set(v, i);

        for (auto w : graph.neighborRange(v)) {
            if (!queue.isProcessed(w)) {
                queue.moveToNewSet(w);
            }
        }
        queue.removeEmptySets();
    }

    return ordering;
}


bool LexBFSChordalGraphRecognition::isZeroFillIn(const PerfectEliminationOrdering& ordering) const {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();

    std::vector<std::vector<NetworKit::node>> A(bound);
    NetworKit::count original_edges = 0;

    for (auto v : graph.nodeRange()) {
        for (auto w : graph.neighborRange(v)) {
            if (ordering.alpha[v] < ordering.alpha[w]) {
                A[v].push_back(w);
                original_edges++;
            }
        }
    }

    std::vector<uint8_t> test(n + 1, false);
    NetworKit::count g_star_edges = 0;

    for (NetworKit::index i = 1; i < n; i++) {
        NetworKit::index k = n;
        NetworKit::node v = ordering.alpha_inv[i];

        std::vector<NetworKit::node> unique_A;

        for (NetworKit::node w : A[v]) {
            NetworKit::index alpha_w = ordering.alpha[w];

            if (!test[alpha_w]) {
                test[alpha_w] = true;
                unique_A.push_back(w);

                k = std::min(k, alpha_w);
            }
        }

        g_star_edges += unique_A.size();

        if (g_star_edges > original_edges) return false;

        NetworKit::node m_v = ordering.alpha_inv[k];

        for (NetworKit::node w : unique_A) {
            test[ordering.alpha[w]] = false;
            if (w != m_v) A[m_v].push_back(w);
        }
    }

    return g_star_edges == original_edges;
}

void LexBFSChordalGraphRecognition::run() {
    hasRun = true;

    if (graph.numberOfNodes() <= 1) {
        is_chordal = State::CHORDAL;
        if (graph.numberOfNodes() == 1) {
            peo = PerfectEliminationOrdering(1, graph.upperNodeIdBound());
            peo.set(*graph.nodeRange().begin(), 1);
        }
        return;
    }

    peo = getLexBFSOrdering();
    is_chordal = isZeroFillIn(peo) ? State::CHORDAL : State::NOT_CHORDAL;
}

} /* namespace Koala */
