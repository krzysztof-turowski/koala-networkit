#include<ranges>

#include <dominatingset/FominKratschWoeginger.hpp>
#include <dominatingset/SchiermeyerMDS.hpp>
#include <tuple>

class RecursiveFKW {
    std::set<NetworKit::node> &free;
    std::set<NetworKit::node> &bounded;
    std::set<NetworKit::node> &required;
    std::set<NetworKit::node> &degreeOne;
    std::set<NetworKit::node> &degreeTwo;
    std::vector<std::set<NetworKit::node>> &neighborhood;

public:
    RecursiveFKW(
        std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bounded,
        std::set<NetworKit::node> &required,
        std::set<NetworKit::node> &degreeOne,
        std::set<NetworKit::node> &degreeTwo,
        std::vector<std::set<NetworKit::node>> &neighborhood)
        : free(free), bounded(bounded), required(required), degreeOne(degreeOne), degreeTwo(degreeTwo), neighborhood(neighborhood) {}
    std::vector<bool> run();
    void onDegreeDecrement(NetworKit::node updated);
    void onDegreeIncrement(NetworKit::node updated);
    void forgetVertexInNeighbors(NetworKit::node vertex);
    void retrieveVertexInNeighbors(NetworKit::node vertex);
    bool forgetVertex(NetworKit::node forced);
    void retrieveVertex(NetworKit::node vertex, bool isFree);
    std::tuple<std::vector<NetworKit::node>, bool> addToTheSolution(NetworKit::node forced);
    void removeFromTheSolution(NetworKit::node vertex, std::vector<NetworKit::node> &movedVertices, bool isVertexFree);
};

FominKratschWoegingerMDS::FominKratschWoegingerMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

void FominKratschWoegingerMDS::run() {
    std::set<NetworKit::node> free;
    std::set<NetworKit::node> bounded;
    std::set<NetworKit::node> required;
    std::set<NetworKit::node> degreeOne;
    std::set<NetworKit::node> degreeTwo;
    std::vector<std::set<NetworKit::node>> neighborhood(G->numberOfNodes());

    G->forNodes([&bounded, &neighborhood, &degreeOne, &degreeTwo, this](NetworKit::node u) {
        bounded.insert(u);
        G->forNeighborsOf(u, [u, &neighborhood](NetworKit::node neighbor) {neighborhood[u].insert(neighbor);});
        if (G->degree(u) == 1) degreeOne.insert(u);
        else if (G->degree(u) == 2) degreeTwo.insert(u);
    });
    dominatingSet = RecursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood).run();
    hasRun = true;
}

std::vector<bool> findMODSWhenDegreeAtLeast3(const NetworKit::Graph &G, const std::set<NetworKit::node> &free, const std::set<NetworKit::node> &bounded) {
    std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);
    for (size_t i = 1; 8 * i <= 3 * (free.size() + bounded.size()); i++) {
        std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);

        auto [found, choices] = SizedChoiceSearcher([&G, &bounded](const std::vector<NetworKit::node> &arg) {return isOptionalDominatingSet(G, arg, bounded);}, possibilities, i).search();
        std::vector<bool> solution(G.numberOfNodes());
        if (found) {
            for (auto u : choices) {
                solution.at(u) = true;
            }
            return solution;
        }
    }
    throw std::invalid_argument("findMODSWhenDegreeAtLeast3 arguments are corrupted");
}

void RecursiveFKW::onDegreeDecrement(NetworKit::node updated) {
    auto currentDegree = neighborhood.at(updated).size();
    if (currentDegree == 0) {
        degreeOne.erase(updated);
    } else if (currentDegree == 1) {
        degreeOne.insert(updated);
        degreeTwo.erase(updated);
    } else if (currentDegree == 2) {
        degreeTwo.insert(updated);
    }
}

void RecursiveFKW::onDegreeIncrement(NetworKit::node updated) {
    auto currentDegree = neighborhood.at(updated).size();
    if (currentDegree == 1) {
        degreeOne.insert(updated);
    } else if (currentDegree == 2) {
        degreeOne.erase(updated);
        degreeTwo.insert(updated);
    } else if (currentDegree == 3) {
        degreeTwo.erase(updated);
    }
}

void RecursiveFKW::forgetVertexInNeighbors(NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).erase(vertex);
        onDegreeDecrement(neighbor);
    }
}
void RecursiveFKW::retrieveVertexInNeighbors(NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).insert(vertex);
        onDegreeIncrement(neighbor);
    }
}

bool RecursiveFKW::forgetVertex(NetworKit::node vertex) {
    bool isFree = free.contains(vertex);
    if (isFree) {
        free.erase(vertex);
    } else {
        bounded.erase(vertex);
    }
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degreeOne.erase(vertex);
    } else if (degree == 2) {
        degreeTwo.erase(vertex);
    }
    forgetVertexInNeighbors(vertex);
    return isFree;
}

void RecursiveFKW::retrieveVertex(NetworKit::node vertex, bool isFree) {
    retrieveVertexInNeighbors(vertex);
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degreeOne.insert(vertex);
    } else if (degree == 2) {
        degreeTwo.insert(vertex);
    }
    if (isFree) {
        free.insert(vertex);
    } else {
        bounded.insert(vertex);
    }
}

std::tuple<std::vector<NetworKit::node>, bool> RecursiveFKW::addToTheSolution(NetworKit::node forced) {
    std::vector<NetworKit::node> movedVertices{};
    for (auto neighbor : neighborhood.at(forced)) {
        if (bounded.contains(neighbor)) {
            bounded.erase(neighbor);
            free.insert(neighbor);
            movedVertices.push_back(neighbor);
        }
    }
    bool isFree = forgetVertex(forced);

    required.insert(forced);
    return {movedVertices, isFree};
}

void RecursiveFKW::removeFromTheSolution(NetworKit::node vertex, std::vector<NetworKit::node> &movedVertices, bool isVertexFree) {
    required.erase(vertex);
    retrieveVertex(vertex, isVertexFree);
    for (auto neighbor : std::views::reverse(movedVertices)) {
        free.erase(neighbor);
        bounded.insert(neighbor);
    }
}

std::vector<bool> RecursiveFKW::run() {
    if (!degreeOne.empty()) {
        NetworKit::node u = *degreeOne.begin();
        NetworKit::node unique = *neighborhood.at(u).begin();
        std::vector<bool> solution;

        bool isFree = forgetVertex(u);

        if (isFree) {
            solution = run();
        } else {
            auto status = addToTheSolution(unique);
            solution = run();
            removeFromTheSolution(unique, std::get<0>(status), std::get<1>(status));
        }
        retrieveVertex(u, isFree);
        return solution;
    }
    if (!degreeTwo.empty()) {
        NetworKit::node v = *degreeTwo.begin();
        std::set<NetworKit::node>::iterator it = neighborhood.at(v).begin();
        NetworKit::node u1 = *(it);
        NetworKit::node u2 = *(++it);

        bool isFree = free.contains(v);
        
        bool is1 = forgetVertex(v);
        auto [moved1, isFree1] = addToTheSolution(u1);
        auto solution1 = run();

        removeFromTheSolution(u1, moved1, isFree1);
        retrieveVertex(v, is1);

        bool is2 = forgetVertex(u1);
        bool is3 = forgetVertex(u2);
        auto [moved2, isFree2] = addToTheSolution(v);

        auto solution2 = run();

        removeFromTheSolution(v, moved2, isFree2);
        retrieveVertex(u2, is3);
        retrieveVertex(u1, is2);

        std::vector<bool> solution3;
        bool is4 = forgetVertex(v);
        if (isFree) {
            solution3 = run();
        } else {
            auto [moved3, isFree3] = addToTheSolution(u2);
            solution3 = run();
            removeFromTheSolution(u2, moved3, isFree3);
        }
        retrieveVertex(v, is4);
        return smallerCardinalitySet(smallerCardinalitySet(solution1, solution2), solution3);
    }
    std::set<NetworKit::node> boundedIsolatedVertices;
    for (int i = 0; i < neighborhood.size(); i++) {
        if (neighborhood.at(i).empty() && bounded.contains(i)) {
            bounded.erase(i);
            boundedIsolatedVertices.insert(i);
        }
    }
    
    std::vector<bool> solution;
    if (bounded.empty()) {
        solution = std::vector<bool>(neighborhood.size());
    } else {
        NetworKit::Graph graph(neighborhood.size());
        for (int i = 0; i < neighborhood.size(); i++) {
            if (required.contains(i)) continue;
            for (auto e : neighborhood.at(i)) {
                if (i < e) {
                    graph.addEdge(i, e);
                }
            }
        }
        solution = findMODSWhenDegreeAtLeast3(graph, free, bounded);
    }
    for (auto u : required) {
        solution.at(u) = true;
    }
    for (auto u : boundedIsolatedVertices) {
        solution.at(u) = true;
    }
    for (auto u : boundedIsolatedVertices) {
        bounded.insert(u);
    }
    return solution;
}
