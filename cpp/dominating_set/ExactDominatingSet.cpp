/*
 * ExactDominatingSet.hpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#include <cassert>
#include <map>
#include <ranges>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <dominating_set/ExactDominatingSet.hpp>

namespace Koala {

bool is_optional_dominating_set(
        const NetworKit::Graph &G,
        const std::vector<NetworKit::node> &choices,
        std::set<NetworKit::node> bounded) {
    for (auto u : choices) {
        bounded.erase(u);
        G.forNeighborsOf(u, [&bounded](NetworKit::node v) { bounded.erase(v); });
    }
    return bounded.empty();
}

std::vector<NetworKit::node> merge(
        const std::set<NetworKit::node> &A, const std::set<NetworKit::node> &B) {
    std::vector<NetworKit::node> merged(A.size() + B.size());
    std::copy(A.begin(), A.end(), merged.begin());
    std::copy(B.begin(), B.end(), merged.begin() + A.size());
    return merged;
}

///////////////
class SizedChoiceSearcher {
    const std::function<bool(const std::vector<NetworKit::node>)> &verifier;
    const std::vector<NetworKit::node> &possibilities;
    int size;
    bool recursive(std::vector<NetworKit::node> &choices, NetworKit::node decideOn, int left);
 public:
    SizedChoiceSearcher(
        const std::function<bool(const std::vector<NetworKit::node>)> &verifier,
        const std::vector<NetworKit::node> &possibilities,
        int size)
        : verifier(verifier), possibilities(possibilities), size(size) {}
    std::tuple<bool, std::vector<NetworKit::node>> search();
};

std::tuple<bool, std::vector<NetworKit::node>> SizedChoiceSearcher::search() {
    std::vector<NetworKit::node> choices{};
    return {recursive(choices, 0, size), choices};
}

bool SizedChoiceSearcher::recursive(
        std::vector<NetworKit::node> &choices,
        NetworKit::node decideOn,
        int left) {
    if (decideOn == possibilities.size()) {
        return verifier(choices);
    }
    if (left > 0) {
        choices.push_back(possibilities.at(decideOn));
        if (recursive(choices, decideOn + 1, left - 1)) {
            return true;
        }
        choices.pop_back();
    }
    return (decideOn + left < possibilities.size() && recursive(choices, decideOn + 1, left));
}

class BigMODSSolver {
    const NetworKit::Graph &graph;
    const std::function<std::vector<bool>(const std::set<NetworKit::node>)>& evaluator;
    std::vector<NetworKit::node> &possibilities;
    std::set<NetworKit::node> choices{};
    std::set<NetworKit::node> closedNeighborhood{};
    void recurse(int decideOn, std::vector<bool> &bestSolution);
 public:
    BigMODSSolver(
        const NetworKit::Graph &graph,
        const std::function<std::vector<bool>(const std::set<NetworKit::node>)>& evaluator,
        std::vector<NetworKit::node> &possibilities)
        : graph(graph), evaluator(evaluator), possibilities(possibilities) {}
    std::vector<bool> run();
};

std::vector<bool> BigMODSSolver::run() {
    std::vector<bool> bestSolution(graph.numberOfNodes());
    for (auto e : possibilities) {
        bestSolution.at(e) = true;
    }
    recurse(0, bestSolution);
    return bestSolution;
}

void BigMODSSolver::recurse(int decideOn, std::vector<bool> &bestSolution) {
    if (decideOn == possibilities.size()) {
        if (closedNeighborhood.size() < 3 * choices.size()
            || closedNeighborhood.size() >= 3 * (choices.size() + 1)) {
            return;
        }
        for (auto e : possibilities) {
            if (choices.contains(e)) {
                continue;
            }
            std::set<int> addedNeighbors;
            if (!closedNeighborhood.contains(e)) {
                addedNeighbors.insert(e);
            }
            graph.forNeighborsOf(e, [&addedNeighbors, this](NetworKit::node neighbor) {
                if (!closedNeighborhood.contains(neighbor)) {
                    addedNeighbors.insert(neighbor);
                }});
            if (closedNeighborhood.size() + addedNeighbors.size() >= 3 * (choices.size() + 1)) {
                return;
            }
        }
        auto optionalDominatingSet = evaluator(choices);
        int optionalDominatingSetSize = std::count(optionalDominatingSet.begin(), optionalDominatingSet.end(), true);
        int bestSolutionSize = std::count(bestSolution.begin(), bestSolution.end(), true);
        if (optionalDominatingSetSize < bestSolutionSize) {
            bestSolution = optionalDominatingSet;
        }
        return;
    }
    recurse(decideOn + 1, bestSolution);

    if (3 * (choices.size() + 1) > possibilities.size()) {
        return;
    }
    std::set<int> addedNeighbors;
    NetworKit::node next = possibilities.at(decideOn);
    if (!closedNeighborhood.contains(next)) {
        addedNeighbors.insert(next);
    }
    graph.forNeighborsOf(
        next,
        [&addedNeighbors, this] (NetworKit::node neighbor) {
            if (!closedNeighborhood.contains(neighbor)) {
                addedNeighbors.insert(neighbor);
            }
        });

    choices.insert(next);
    for (auto e : addedNeighbors) {
        closedNeighborhood.insert(e);
    }
    recurse(decideOn + 1, bestSolution);
    choices.erase(next);
    for (auto e : addedNeighbors) {
        closedNeighborhood.erase(e);
    }
}

NetworKit::Graph buildFromVectorSetRepresentation(
        const std::vector<std::set<NetworKit::node>> &neighbors) {
    NetworKit::Graph graph(neighbors.size());
    for (NetworKit::node i = 0; i < neighbors.size(); i++) {
        for (auto e : neighbors.at(i)) {
            assert(neighbors.at(e).contains(i));
            if (i < e) {
                graph.addEdge(i, e);
            }
        }
    }
    return graph;
}

NetworKit::Graph core(
        const NetworKit::Graph &G,
        std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bounded,
        std::set<NetworKit::node> &required) {
    std::vector<std::set<NetworKit::node>> intermediate(G.numberOfNodes());
    G.forEdges([&intermediate](NetworKit::node u, NetworKit::node v) {
        intermediate[u].insert(v);
        intermediate[v].insert(u);
    });
    bool process = true;
    while (process) {
        process = false;
        for (auto e : bounded) {
            if (intermediate[e].empty()) {
                process = true;
                required.insert(e);
            }
        }
        std::erase_if(
            bounded,
            [&intermediate](NetworKit::node u) {
                return intermediate[u].empty();
            });
        std::set<NetworKit::node> rule2change;
        for (auto e : bounded) {
            for (auto nei : intermediate[e]) {
                if (required.contains(nei)) {
                    process = true;
                    rule2change.insert(e);
                    break;
                }
            }
        }
        for (auto e : rule2change) {
            free.insert(e);
        }
        std::erase_if(bounded, [&rule2change](NetworKit::node u) {return rule2change.contains(u);});
        for (int i = 0; i < intermediate.size(); i++) {
            if (free.contains(i)) {
                auto numberOfErased = std::erase_if(
                    intermediate[i],
                    [&free, &required](NetworKit::node e) {
                        return free.contains(e) || required.contains(e);
                    });
                if (numberOfErased > 0) {
                    process = true;
                }
            } else if (required.contains(i)) {
                if (!intermediate[i].empty()) {
                    process = true;
                    for (auto e : intermediate[i]) {
                        if (bounded.contains(e)) {
                            intermediate.at(e).erase(i);
                        }
                    }
                    intermediate[i].clear();
                }
            }
        }
        std::set<NetworKit::node> rule4change;
        for (auto e : free) {
            int boundCount = 0;
            for (auto nei : intermediate[e]) {
                if (bounded.contains(nei)) {
                    boundCount++;
                }
                if (boundCount > 1) {
                    break;
                }
            }
            if (boundCount <= 1) {
                process = true;
                rule4change.insert(e);
            }
        }
        for (auto e : rule4change) {
            for (auto nei : intermediate[e]) {
                intermediate[nei].erase(e);
            }
            intermediate[e].clear();
        }
        std::erase_if(free, [&rule4change](NetworKit::node u) {return rule4change.contains(u);});
        for (auto e : bounded) {
            if (intermediate[e].size() == 1 && !required.contains(e)) {
                process = true;
                auto unique = *intermediate[e].begin();
                free.erase(unique);
                required.insert(unique);
            }
        }
        std::erase_if(bounded, [&required](NetworKit::node u) {return required.contains(u);});
    }
    return buildFromVectorSetRepresentation(intermediate);
}

std::vector<bool> matchingMODS(
        const NetworKit::Graph &graph,
        const std::set<NetworKit::node> &S,
        const std::set<NetworKit::node> &free,
        const std::set<NetworKit::node> &bounded) {
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
    auto SPrim(S);
    auto TPrim(free);
    auto UPrim(bounded);
    auto coreGraph = core(graph, TPrim, UPrim, SPrim);
    boost_graph_t H(graph.numberOfNodes());
    std::map<std::tuple<NetworKit::node, NetworKit::node>, NetworKit::node> owners;
    for (auto u : UPrim) {
        NetworKit::node v = coreGraph.getIthNeighbor(u, 0);
        if (v != NetworKit::none && u < v) {
            boost::add_edge(u, v, H);
            owners.emplace(std::tuple<NetworKit::node, NetworKit::node>(u, v), u);
        }
    }

    for (auto u : TPrim) {
        auto v1 = coreGraph.getIthNeighbor(u, 0), v2 = coreGraph.getIthNeighbor(u, 1);
        assert(v1 != NetworKit::none);
        assert(v2 != NetworKit::none);
        if (v1 > v2) {
            std::swap(v1, v2);
        }
        if (!owners.contains(std::tuple<NetworKit::node, NetworKit::node>(v1, v2))) {
            boost::add_edge(v1, v2, H);
            owners.emplace(std::tuple<NetworKit::node, NetworKit::node>(v1, v2), u);
        }
    }
    std::vector<boost::graph_traits<boost_graph_t>::vertex_descriptor> mate(graph.numberOfNodes());
    assert(boost::checked_edmonds_maximum_cardinality_matching(H, &mate[0]));

    std::vector<bool> optionalDominatingSet(graph.numberOfNodes());
    for (auto e : UPrim) {
        if (mate[e] == NetworKit::none) {
            optionalDominatingSet[e] = true;
        } else if (e < mate[e]) {
            auto owner = owners.at(std::tuple<NetworKit::node, NetworKit::node>(e, mate[e]));
            optionalDominatingSet[owner] = true;
        }
    }
    for (auto e : SPrim) {
        optionalDominatingSet[e] = true;
    }
    return optionalDominatingSet;
}

///////////////

void SchiermeyerDominatingSet::run() {
    hasRun = true;
    dominating_set.resize(graph->numberOfNodes());
    graph->forNodes([this](NetworKit::node u) {
        bounded.insert(u);
    });
    NetworKit::Graph core_graph = core(*graph, free, bounded, required);
    if (bounded.empty()) {
        graph->forNodes([this](NetworKit::node u) {
            dominating_set.at(u) = required.contains(u);
        });
    }
    if (find_small_MODS(core_graph)) {
        return;
    }
    find_big_MODS(core_graph);
}

bool SchiermeyerDominatingSet::find_small_MODS(const NetworKit::Graph &G) {
    auto possibilities = merge(free, bounded);
    for (int i = 1; 3 * i <= possibilities.size(); i++) {
        auto [found, choices] = SizedChoiceSearcher(
                [&G, this](const std::vector<NetworKit::node> &arg) {
                    return is_optional_dominating_set(G, arg, bounded);
                }, possibilities, i).search();
        if (found) {
            G.forNodes([this](NetworKit::node u) {
                dominating_set.at(u) = required.contains(u);
            });
            for (auto u : choices) {
                dominating_set.at(u) = true;
            }
            return true;
        }
    }
    return false;
}

void SchiermeyerDominatingSet::find_big_MODS(const NetworKit::Graph &G) {
    auto possibilities = merge(free, bounded);
    dominating_set = BigMODSSolver(
        G,
        [&G, this](const std::set<NetworKit::node> &S) {
            return matchingMODS(G, S, free, bounded);
        }, possibilities).run();
    for (auto e : required) {
        dominating_set.at(e) = true;
    }
}

//////////////////////

std::vector<bool> &smaller(std::vector<bool> &lhs, std::vector<bool> &rhs) {
    if (std::count(lhs.begin(), lhs.end(), true) < std::count(rhs.begin(), rhs.end(), true)) {
        return lhs;
    } else {
        return rhs;
    }
}

void FominKratschWoegingerDominatingSet::run() {
    hasRun = true;
    neighborhood.resize(graph->numberOfNodes());
    graph->forNodes([this](NetworKit::node u) {
        bounded.insert(u);
        graph->forNeighborsOf(u, [u, this](NetworKit::node neighbor) {
            neighborhood[u].insert(neighbor);
        });
        if (graph->degree(u) == 1) {
            degreeOne.insert(u);
        } else if (graph->degree(u) == 2) {
            degreeTwo.insert(u);
        }
    });
    dominating_set = recurse();
}

std::vector<bool> findMODSWhenDegreeAtLeast3(
        const NetworKit::Graph &G,
        const std::set<NetworKit::node> &free,
        const std::set<NetworKit::node> &bounded) {
    std::vector<NetworKit::node> possibilities = merge(free, bounded);
    for (size_t i = 1; 8 * i <= 3 * possibilities.size(); i++) {
        auto [found, choices] = SizedChoiceSearcher(
            [&G, &bounded](const std::vector<NetworKit::node> &arg) {
                return is_optional_dominating_set(G, arg, bounded);
            }, possibilities, i).search();
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

void FominKratschWoegingerDominatingSet::onDegreeDecrement(NetworKit::node updated) {
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

void FominKratschWoegingerDominatingSet::onDegreeIncrement(NetworKit::node updated) {
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

void FominKratschWoegingerDominatingSet::forgetVertexInNeighbors(NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).erase(vertex);
        onDegreeDecrement(neighbor);
    }
}
void FominKratschWoegingerDominatingSet::retrieveVertexInNeighbors(NetworKit::node vertex) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).insert(vertex);
        onDegreeIncrement(neighbor);
    }
}

bool FominKratschWoegingerDominatingSet::forgetVertex(NetworKit::node vertex) {
    bool is_free = free.contains(vertex);
    if (is_free) {
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
    return is_free;
}

void FominKratschWoegingerDominatingSet::retrieveVertex(NetworKit::node vertex, bool is_free) {
    retrieveVertexInNeighbors(vertex);
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degreeOne.insert(vertex);
    } else if (degree == 2) {
        degreeTwo.insert(vertex);
    }
    if (is_free) {
        free.insert(vertex);
    } else {
        bounded.insert(vertex);
    }
}

std::tuple<std::vector<NetworKit::node>, bool> FominKratschWoegingerDominatingSet::addToTheSolution(
        NetworKit::node forced) {
    std::vector<NetworKit::node> movedVertices{};
    for (auto neighbor : neighborhood.at(forced)) {
        if (bounded.contains(neighbor)) {
            bounded.erase(neighbor);
            free.insert(neighbor);
            movedVertices.push_back(neighbor);
        }
    }
    bool is_free = forgetVertex(forced);
    required.insert(forced);
    return {movedVertices, is_free};
}

void FominKratschWoegingerDominatingSet::removeFromTheSolution(
        NetworKit::node vertex, std::vector<NetworKit::node> &movedVertices, bool isVertexFree) {
    required.erase(vertex);
    retrieveVertex(vertex, isVertexFree);
    for (auto neighbor : std::views::reverse(movedVertices)) {
        free.erase(neighbor);
        bounded.insert(neighbor);
    }
}

std::vector<bool> FominKratschWoegingerDominatingSet::recurse() {
    if (!degreeOne.empty()) {
        NetworKit::node u = *degreeOne.begin(), unique = *neighborhood.at(u).begin();
        std::vector<bool> solution;
        bool is_free = forgetVertex(u);
        if (is_free) {
            solution = recurse();
        } else {
            auto [moved1, is_free1] = addToTheSolution(unique);
            solution = recurse();
            removeFromTheSolution(unique, moved1, is_free1);
        }
        retrieveVertex(u, is_free);
        return solution;
    }
    if (!degreeTwo.empty()) {
        NetworKit::node v = *degreeTwo.begin();
        std::set<NetworKit::node>::iterator it = neighborhood.at(v).begin();
        NetworKit::node u1 = *it, u2 = *(++it);

        bool is_free = free.contains(v), is1 = forgetVertex(v);
        auto [moved1, is_free1] = addToTheSolution(u1);
        auto solution1 = recurse();
        removeFromTheSolution(u1, moved1, is_free1);
        retrieveVertex(v, is1);

        bool is2 = forgetVertex(u1), is3 = forgetVertex(u2);
        auto [moved2, is_free2] = addToTheSolution(v);
        auto solution2 = recurse();
        removeFromTheSolution(v, moved2, is_free2);
        retrieveVertex(u2, is3), retrieveVertex(u1, is2);

        std::vector<bool> solution3;
        bool is4 = forgetVertex(v);
        if (is_free) {
            solution3 = recurse();
        } else {
            auto [moved3, is_free3] = addToTheSolution(u2);
            solution3 = recurse();
            removeFromTheSolution(u2, moved3, is_free3);
        }
        retrieveVertex(v, is4);
        return smaller(smaller(solution1, solution2), solution3);
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
        for (NetworKit::node i = 0; i < neighborhood.size(); i++) {
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

}  /* namespace Koala */
