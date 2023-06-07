#include <coloring/InclusionExclusionVertexColoring.hpp>
#include <iostream>
#include <optional>
#include <vector>
namespace Koala {
IndependentSetChecker::IndependentSetChecker(const NetworKit::Graph& graph)
: graph(std::make_optional(graph)) {
}
bool IndependentSetChecker::isIndependentSet(std::vector<NetworKit::node>& nodes) {
    if (!graph.has_value()) {
        throw std::runtime_error("Graph is not set.");
    }
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
        for (auto j = i + 1; j != nodes.end(); ++j) {
            if (graph->hasEdge(*i, *j)) {
                return false;
            }
        }
    }
    return true;
}
long long int IndependentSetChecker::numberOfIndependentSetsNotIntersectingWith(
std::vector<NetworKit::node>& nodes) {

    if (!graph.has_value()) {
        throw std::runtime_error("Graph is not set.");
    }
    std::unordered_set<int> nodesInX(nodes.begin(), nodes.end());
    std::vector<NetworKit::node> nodesWithoutX;
    for (int i = 0; i < graph->numberOfNodes(); i++) {
        if (nodesInX.find(i) == nodesInX.end())
            nodesWithoutX.push_back(i);
    }
    long long int result = 0;
    for (int i = 0; i < (1 << nodesWithoutX.size()); i++) {
        std::vector<NetworKit::node> subset;
        for (int j = 0; j < nodesWithoutX.size(); j++) {
            if (i & (1 << j)) {
                subset.push_back(nodesWithoutX[j]);
            }
        }
        if (isIndependentSet(subset)) {
            result++;
        }
    }
    return result;
}

InclusionExclusionVertexColoring::InclusionExclusionVertexColoring(const NetworKit::Graph& graph)
: graph(std::make_optional(graph)) {
}

const std::map<NetworKit::node, int> InclusionExclusionVertexColoring::getColoring() const {
    assureFinished();
    std::map<NetworKit::node, int> coloring;
    for (int i = 0; i < graph->numberOfNodes(); i++) {
        coloring[i] = best_solution[i];
    }
    return coloring;
}

bool InclusionExclusionVertexColoring::is_k_colorable(int k) {
    long long int numberOfKCovers = 0;
    IndependentSetChecker independentSetChecker(*graph);
    std::vector<NetworKit::node> X;
    for (int i = 0; i < (1 << graph->numberOfNodes()); i++) {
        X.clear();
        for (int j = 0; j < graph->numberOfNodes(); j++) {
            if (i & (1 << j)) {
                X.push_back(j);
            }
        }
        long long int result = ((X.size() % 2) ? -1 : 1) *
        pow(independentSetChecker.numberOfIndependentSetsNotIntersectingWith(X), k);
        numberOfKCovers += result;
    }
    return numberOfKCovers > 0;
}

void InclusionExclusionVertexColoring::run() {
    int n = graph->numberOfNodes();
    int l = 1;
    int r = n - 1;
    best_solution = std::vector<int>(n, 0);
    if (is_k_colorable(l)) {
        best_solution = std::vector<int>(n, 1);
        chromatic_number = 1;
    } else if (!is_k_colorable(r)) {
        for (int i = 0; i < n; i++) {
            best_solution[i] = i + 1;
        }
        chromatic_number = n;
    } else {
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (is_k_colorable(m)) {
                r = m;
            } else {
                l = m + 1;
            }
        }
        if (is_k_colorable(l)) {
            best_solution = std::vector<int>(n, l);
            chromatic_number = l;
        } else {
            best_solution = std::vector<int>(n, r);
            chromatic_number = r;
        }
    }

    hasRun = true;
}

const int InclusionExclusionVertexColoring::getChromaticNumber() const {
    assureFinished();
    return chromatic_number;
}

} // namespace Koala
