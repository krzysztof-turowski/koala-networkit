/*
 * MaximumCardinalitySearchChordalGraphRecognition.cpp
 *
 *  Created on: 2026-04-14
 *      Author: Mateusz Przebieracz
 */

#include "recognition/ChordalGraphRecognition.hpp"

namespace Koala {


MaximumCardinalitySearchChordalGraphRecognition::SetArray::SetArray(const NetworKit::Graph& graph) {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();

    head.assign(n, NetworKit::none);
    vertices.resize(bound);

    for (auto v : graph.nodeRange()) add(v, INITIAL_SET);
}

void MaximumCardinalitySearchChordalGraphRecognition::SetArray::add(
    NetworKit::node v, NetworKit::count set) {
    vertices[v].next = head[set];
    vertices[v].prev = NetworKit::none;
    if (head[set] != NetworKit::none) {
        vertices[head[set]].prev = v;
    }
    head[set] = v;
}

void MaximumCardinalitySearchChordalGraphRecognition::SetArray::remove(
    NetworKit::node v, NetworKit::count set) {
    if (vertices[v].prev != NetworKit::none) {
        vertices[vertices[v].prev].next = vertices[v].next;
    }
    else {
        head[set] = vertices[v].next;
    }

    if (vertices[v].next != NetworKit::none) {
        vertices[vertices[v].next].prev = vertices[v].prev;
    }

    vertices[v].prev = NetworKit::none;
    vertices[v].next = NetworKit::none;
}

bool MaximumCardinalitySearchChordalGraphRecognition::SetArray::isEmpty(NetworKit::count set) const {
    return head[set] == NetworKit::none;
}

NetworKit::node MaximumCardinalitySearchChordalGraphRecognition::SetArray::getHead(NetworKit::count set) const {
    return head[set];
}


PerfectEliminationOrdering MaximumCardinalitySearchChordalGraphRecognition::getMCSOrdering() const {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();
    constexpr NetworKit::count PROCESSED_SIZE = std::numeric_limits<NetworKit::count>::max();

    PerfectEliminationOrdering ordering(n, bound);

    std::vector<NetworKit::count> size(bound, 0);
    SetArray sets(graph);

    NetworKit::count max_nonempty_set = SetArray::INITIAL_SET;

    for (NetworKit::count i = n; i >= 1; i--) {
        while (max_nonempty_set > 0 && sets.isEmpty(max_nonempty_set)) {
            max_nonempty_set--;
        }

        NetworKit::node v = sets.getHead(max_nonempty_set);
        sets.remove(v, max_nonempty_set);

        ordering.set(v, i);
        size[v] = PROCESSED_SIZE;

        for (auto w : graph.neighborRange(v)) {
            if (size[w] != PROCESSED_SIZE) {
                sets.remove(w, size[w]);
                size[w]++;
                sets.add(w, size[w]);
            }
        }

        max_nonempty_set++;
    }

    return ordering;
}

bool MaximumCardinalitySearchChordalGraphRecognition::isZeroFillIn(const PerfectEliminationOrdering& ordering) const {
    NetworKit::count n = graph.numberOfNodes();
    NetworKit::count bound = graph.upperNodeIdBound();

    std::vector<NetworKit::node> follower(bound);
    std::vector<NetworKit::count> index(bound);

    for (NetworKit::count i = 1; i <= n; i++) {
        NetworKit::node w = ordering.alpha_inv[i];
        follower[w] = w;
        index[w] = i;

        for (auto v : graph.neighborRange(w)) {
            if (ordering.alpha[v] < i) {
                index[v] = i;
                if (follower[v] == v) {
                    follower[v] = w;
                }
            }
        }

        for (auto v : graph.neighborRange(w)) {
            if (ordering.alpha[v] < i && index[follower[v]] < i) {
                return false;
            }
        }
    }

    return true;
}

void MaximumCardinalitySearchChordalGraphRecognition::run() {
    hasRun = true;

    if (graph.numberOfNodes() <= 1) {
        is_chordal = State::CHORDAL;
        return;
    }

    PerfectEliminationOrdering ordering = getMCSOrdering();

    is_chordal = isZeroFillIn(ordering) ? State::CHORDAL : State::NOT_CHORDAL;
}

} /* namespace Koala */
