#include <cassert>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <dominating_set/ExactDominatingSets.hpp>

namespace Koala {

NetworKit::Graph core(
    const NetworKit::Graph &G,
    std::set<NetworKit::node> &free,
    std::set<NetworKit::node> &bounded,
    std::set<NetworKit::node> &required);
std::tuple<bool, std::vector<bool>> findSmallMODS(
    const NetworKit::Graph &G,
    const std::set<NetworKit::node> &free,
    const std::set<NetworKit::node> &bounded,
    const std::set<NetworKit::node> &required);
std::vector<bool> findBigMODS(
    const NetworKit::Graph &G,
    const std::set<NetworKit::node> &free,
    const std::set<NetworKit::node> &bounded,
    const std::set<NetworKit::node> &required);

SchiermeyerMDS::SchiermeyerMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

void SchiermeyerMDS::run() {
    std::set<NetworKit::node> free;
    std::set<NetworKit::node> bounded;
    std::set<NetworKit::node> required;

    G->forNodes([&bounded](NetworKit::node u) {
        bounded.insert(u);
    });
    NetworKit::Graph coreGraph = core(*G, free, bounded, required);
    if (bounded.empty()) {
        dominatingSet = std::vector<bool>(G->numberOfNodes());
        G->forNodes([&required, this](NetworKit::node u) {
            this->dominatingSet.at(u) = required.contains(u);
        });
    }
    auto [found, smallSolution] = findSmallMODS(coreGraph, free, bounded, required);
    if (found) {
        dominatingSet = smallSolution;
    } else {
        dominatingSet = findBigMODS(coreGraph, free, bounded, required);
    }
    hasRun = true;
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

std::vector<NetworKit::node> joinFreeAndBounded(
        const std::set<NetworKit::node> &free,
        const std::set<NetworKit::node> &bounded) {
    std::vector<NetworKit::node> possibilities;
    for (auto e : free) {
        possibilities.push_back(e);
    }
    for (auto e : bounded) {
        possibilities.push_back(e);
    }
    return possibilities;
}

std::tuple<bool, std::vector<bool>> findSmallMODS(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &free,
        const std::set<NetworKit::node> &bounded,
        const std::set<NetworKit::node> &required) {
    const std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);
    for (int i = 1; 3 * i <= free.size() + bounded.size(); i++) {
        auto [found, choices] = SizedChoiceSearcher(
                [&G, &bounded](const std::vector<NetworKit::node> &arg) {
                    return isOptionalDominatingSet(G, arg, bounded);
                },
                possibilities,
                i)
            .search();
        std::vector<bool> solution(G.numberOfNodes());
        if (found) {
            G.forNodes([&required, &solution](NetworKit::node u) {
                solution.at(u) = required.contains(u);
            });
            for (auto u : choices) {
                solution.at(u) = true;
            }
            return {true, solution};
        }
    }
    return {false, std::vector<bool>{}};
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
        if (MinimumDominatingSet::dominatingSetSize(optionalDominatingSet) <
            MinimumDominatingSet::dominatingSetSize(bestSolution)) {
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
        auto v1 = coreGraph.getIthNeighbor(u, 0);
        auto v2 = coreGraph.getIthNeighbor(u, 1);
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

    bool success = boost::checked_edmonds_maximum_cardinality_matching(H, &mate[0]);
    assert(success);

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

std::vector<bool> findBigMODS(
        const NetworKit::Graph &G,
        const std::set<NetworKit::node> &free,
        const std::set<NetworKit::node> &bounded,
        const std::set<NetworKit::node> &required) {
    std::vector<NetworKit::node> possibilities = joinFreeAndBounded(free, bounded);
    auto bestSolution = BigMODSSolver(
            G,
            [&G, &free, &bounded](const std::set<NetworKit::node> S) {
                return matchingMODS(G, S, free, bounded);
            },
            possibilities)
        .run();
    for (auto e : required) {
        bestSolution.at(e) = true;
    }
    return bestSolution;
}

bool isOptionalDominatingSet(
        const NetworKit::Graph &G,
        const std::vector<NetworKit::node> &choices,
        const std::set<NetworKit::node> &bounded) {
    std::set<NetworKit::node> boundedCopy(bounded);
    for (auto e : choices) {
        boundedCopy.erase(e);
        G.forNeighborsOf(
            e,
            [&boundedCopy](NetworKit::node neighbor) {
                boundedCopy.erase(neighbor);
            });
    }
    return boundedCopy.empty();
}
}  /* namespace Koala */
