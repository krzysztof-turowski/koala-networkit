#include<ranges>

#include <dominatingset/FominKratschWoeginger.hpp>
#include <dominatingset/SchiermeyerMDS.hpp>
#include <tuple>
std::vector<bool> recursiveFKW(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &required, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, int depth);

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
    dominatingSet = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, 0);
    hasRun = true;
}

std::vector<bool> findMODSWhenDegreeAtLeast3(const NetworKit::Graph &G, const std::set<NetworKit::node> &free, const std::set<NetworKit::node> &bounded) {
    std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);
    for (size_t i = 1; 8 * i <= 3 * (free.size() + bounded.size()); i++) {
        std::vector<NetworKit::node> choices{};
        std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);

        bool found = recursiveSizedChoiceSearch([&G, &bounded](const std::vector<NetworKit::node> &arg) {return isOptionalDominatingSet(G, arg, bounded);}, possibilities, choices, 0, i);
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

void onDegreeDecrement(std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, NetworKit::node updated, int currentDegree) {
    if (currentDegree == 0) {
        degreeOne.erase(updated);
    } else if (currentDegree == 1) {
        degreeOne.insert(updated);
        degreeTwo.erase(updated);
    } else if (currentDegree == 2) {
        degreeTwo.insert(updated);
    }
}

void onDegreeIncrement(std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, NetworKit::node updated, int currentDegree) {
    if (currentDegree == 1) {
        degreeOne.insert(updated);
    } else if (currentDegree == 2) {
        degreeOne.erase(updated);
        degreeTwo.insert(updated);
    } else if (currentDegree == 3) {
        degreeTwo.erase(updated);
    }
}

void forgetVertexInNeighbors(std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).erase(vertex);
        onDegreeDecrement(degreeOne, degreeTwo, neighbor, neighborhood.at(neighbor).size());
    }
}
void retrieveVertexInNeighbors(std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).insert(vertex);
        onDegreeIncrement(degreeOne, degreeTwo, neighbor, neighborhood.at(neighbor).size());
    }
}

bool forgetVertex(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node vertex) {
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
    forgetVertexInNeighbors(degreeOne, degreeTwo, neighborhood, vertex);
    return isFree;
}

void retrieveVertex(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node vertex, bool isFree) {
    retrieveVertexInNeighbors(degreeOne, degreeTwo, neighborhood, vertex);
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

std::tuple<std::vector<NetworKit::node>, bool> addToTheSolution(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &required, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node forced) {
    std::vector<NetworKit::node> movedVertices{};
    for (auto neighbor : neighborhood.at(forced)) {
        if (bounded.contains(neighbor)) {
            bounded.erase(neighbor);
            free.insert(neighbor);
            movedVertices.push_back(neighbor);
        }
    }
    bool isFree = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, forced);

    required.insert(forced);
    return {movedVertices, isFree};
}

void removeFromTheSolution(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &required, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, NetworKit::node vertex, std::vector<NetworKit::node> &movedVertices, bool isVertexFree) {
    required.erase(vertex);
    retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, vertex, isVertexFree);
    for (auto neighbor : std::views::reverse(movedVertices)) {
        free.erase(neighbor);
        bounded.insert(neighbor);
    }
}

void log(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &required, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, int depth) {
    // std::cout << "log " << depth << "\n";
    // std::cout << "free ";
    // for (auto e : free) {
    //     std::cout << e << " ";
    // }std::cout << "\t";
    // std::cout << "bounded ";
    // for (auto e : bounded) {
    //     std::cout << e << " ";
    // }std::cout << "\t";
    // std::cout << "required ";
    // for (auto e : required) {
    //     std::cout << e << " ";
    // }std::cout << "\t";
    // std::cout << "one ";
    // for (auto e : degreeOne) {
    //     std::cout << e << " ";
    // }std::cout << "\t";
    // std::cout << "two ";
    // for (auto e : degreeTwo) {
    //     std::cout << e << " ";
    // }std::cout << "\n";
    // std::cout << "neighborhood";
    // for (auto e : neighborhood) {
    //     std::cout << "[";
    //     for (auto f : e) {
    //         std::cout << f << " ";
    //     }std::cout << "]\t";
    // }std::cout << "\n";
}

std::vector<bool> recursiveFKW(std::set<NetworKit::node> &free, std::set<NetworKit::node> &bounded, std::set<NetworKit::node> &required, std::set<NetworKit::node> &degreeOne, std::set<NetworKit::node> &degreeTwo, std::vector<std::set<NetworKit::node>> &neighborhood, int depth) {
    log(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth);
    if (!degreeOne.empty()) {
        NetworKit::node u = *degreeOne.begin();
        NetworKit::node unique = *neighborhood.at(u).begin();
        std::vector<bool> solution;

        bool isFree = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u);

        if (isFree) {
            solution = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);
        } else {
            auto status = addToTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, unique);
            solution = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);
            removeFromTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, unique, std::get<0>(status), std::get<1>(status));
        }
        retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u, isFree);
        log(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth);
        return solution;
    }
    if (!degreeTwo.empty()) {
        NetworKit::node v = *degreeTwo.begin();
        std::set<NetworKit::node>::iterator it = neighborhood.at(v).begin();
        NetworKit::node u1 = *(it);
        NetworKit::node u2 = *(++it);

        bool isFree = free.contains(v);
        
        bool is1 = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, v);
        auto [moved1, isFree1] = addToTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, u1);
        auto solution1 = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);

        removeFromTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, u1, moved1, isFree1);
        retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, v, is1);

        log(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth);

        bool is2 = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u1);
        bool is3 = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u2);
        auto [moved2, isFree2] = addToTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, v);

        auto solution2 = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);

        removeFromTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, v, moved2, isFree2);
        retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u2, is3);
        retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, u1, is2);

        log(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth);

        std::vector<bool> solution3;
        bool is4 = forgetVertex(free, bounded, degreeOne, degreeTwo, neighborhood, v);
        if (isFree) {
            solution3 = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);
        } else {
            auto [moved3, isFree3] = addToTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, u2);
            solution3 = recursiveFKW(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth + 1);
            removeFromTheSolution(free, bounded, required, degreeOne, degreeTwo, neighborhood, u2, moved3, isFree3);
        }
        retrieveVertex(free, bounded, degreeOne, degreeTwo, neighborhood, v, is4);
        log(free, bounded, required, degreeOne, degreeTwo, neighborhood, depth);
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
