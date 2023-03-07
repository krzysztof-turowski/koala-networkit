
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

LightWeightGraph::LightWeightGraph(const NetworKit::Graph &graph)
        : adj(std::vector<std::vector<int>>(graph.numberOfNodes()))
        , hidden(std::vector<bool>(graph.numberOfNodes(), false)) {
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    });
}

void LightWeightGraph::hide(int v) {
    hidden[v] = true;
}

void LightWeightGraph::unhide(int v) {
    hidden[v] = false;
}

void LightWeightGraph::hide(std::vector<int>& v) {
    for (int i = 0; i < v.size(); ++i) {
        hidden[v[i]] = true;
    }
}

void LightWeightGraph::unhide(std::vector<int>& v) {
    for (int i = 0; i < v.size(); ++i) {
        hidden[v[i]] = false;
    }
}

std::vector<int> LightWeightGraph::n(int v) {
    std::vector<int> neighbors;
    neighbors.push_back(v);
    for (int i = 0; i < adj[v].size(); ++i) {
        if (!hidden[adj[v][i]]) {
            neighbors.push_back(adj[v][i]);
        }
    }
    return neighbors;
}

int LightWeightGraph::lowestDegVerticle() {
    int lowestDeg = adj.size();
    int lowestDegIndex;

    for (int i = 0; i < adj.size(); ++i) {
        if (!hidden[i]) {
            int deg = 0;
            for (int j = 0; j < adj[i].size(); ++j) {
                if (!hidden[adj[i][j]]) {
                    ++deg;
                }
            }
            if (deg < lowestDeg) {
                lowestDeg = deg;
                lowestDegIndex = i;
            }
        }
    }

    return lowestDegIndex;
}

bool LightWeightGraph::isEmpty() {
    for (int i = 0; i < hidden.size(); ++i) {
        if (!hidden[i]) {
            return false;
        }
    }
    return true;
}

void LightWeightGraph::print() {
    for (int i = 0; i < adj.size(); ++i) {
        if (!hidden[i]) {
            std::cout << i << ": ";
            for (int j = 0; j < adj[i].size(); ++j) {
                if (!hidden[adj[i][j]]) {
                    std::cout << adj[i][j] << " ";
                }
            }
            std::cout << std::endl;
        }
    }
}


SimpleIndependentSet::SimpleIndependentSet(const NetworKit::Graph &graph) 
    : graph(std::make_optional(graph))
    , lightGraph(graph) { }

const std::map<NetworKit::node, bool>& SimpleIndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void BruteForceIndependentSet::run() {
    if (graph->numberOfNodes() == 0) {
        return;
    }

    std::vector<bool> testSet;
    graph->forNodes([&](NetworKit::node v) {
        testSet.push_back(false);
        independentSet[v] = false;
    });

    int best = 0;
    unsigned long long max = 1;
    max <<= graph->numberOfNodes();

    for (unsigned long long binary = 1; binary != max; ++binary) {
        unsigned long long tmp = binary;
        int testSetSize = 0;
        for (int i = 0; i < testSet.size(); ++i) {
            testSet[i] = (tmp % 2);
            testSetSize += (tmp % 2);
            tmp >>= 1;
        }

        if (testSetSize <= best) {
            continue;
        }

        bool bad = false;
        graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
            if(testSet[u] && testSet[v]) {
                bad = true;
                return;
            }
        });

        if (!bad) {
            for (int i = 0; i < testSet.size(); ++i) {
                independentSet[i] = testSet[i];
                best = testSetSize;
            }
        }
    }
    hasRun = true;
}

static void print_vector(std::vector<int>& v) {
    for (int i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
    }
}

void Mis1IndependentSet::run() {
    std::vector<int> result = recursive();
    graph->forNodes([&](NetworKit::node v) {
        independentSet[v] = false;
    });
    for (int i = 0; i < result.size(); ++i) {
        independentSet[result[i]] = true;
    }
    hasRun = true;
}

std::vector<int> Mis1IndependentSet::recursive() {
    if (lightGraph.isEmpty()) {
        return {};
    }
    int v = lightGraph.lowestDegVerticle();
    std::vector<int> oneMustBeInSet = lightGraph.n(v);

    int chosenToSet;
    std::vector<int> best;
    for (int i = 0; i < oneMustBeInSet.size(); ++i) {
        std::vector<int> neighbors = lightGraph.n(oneMustBeInSet[i]);
        lightGraph.hide(neighbors);
        
        std::vector<int> result = recursive();
        if (result.size() >= best.size()) {
            best = result;
            chosenToSet = oneMustBeInSet[i];
        }

        lightGraph.unhide(neighbors);
    }
    best.push_back(chosenToSet);
    return best;
}


void Mis2IndependentSet::run() {

}

void Mis3IndependentSet::run() {

}

void Mis4IndependentSet::run() {

}

void Mis5IndependentSet::run() {

}

} /* namespace Koala */
