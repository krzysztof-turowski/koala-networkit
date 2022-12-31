
/*
 * SimpleIndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independentSet/SimpleIndependentSet.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace Koala {

SimpleIndependentSet::SimpleIndependentSet(
        const NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::map<NetworKit::node, bool>& SimpleIndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void BruteForceIndependentSet::run() {
    if (graph->numberOfNodes() == 0) {
        return;
    }

    std::vector<std::pair<NetworKit::node, bool>> test;
    graph->forNodes([&](NetworKit::node v) {
        test.push_back({v, false});
        independentSet[v] = false;
    });

    int best = 0;
    unsigned long long max = 1;
    max <<= graph->numberOfNodes();

    for (unsigned long long binary = 1; binary != max; ++binary) {
        unsigned long long tmp = binary;
        int testSize = 0;
        for (int i = 0; i < test.size(); ++i) {
            test[i].second = (tmp % 2);
            testSize += (tmp % 2);
            tmp >>= 1;
        }

        if (testSize <= best) {
            continue;
        }

        bool bad = false;
        graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
            if(test[u].second && test[v].second) {
                bad = true;
                return;
            }
        });

        if (!bad) {
            for (int i = 0; i < test.size(); ++i) {
                independentSet[i] = test[i].second;
                best = testSize;
            }
        }
    }
    hasRun = true;
}

} /* namespace Koala */
