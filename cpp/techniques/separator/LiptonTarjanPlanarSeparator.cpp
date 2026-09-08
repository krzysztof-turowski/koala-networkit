#include "techniques/separator/LiptonTarjanPlanarSeparator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/graph/filtered_graph.hpp>
#include <shortest_path/planar/SuitableRDivision.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>

#include "techniques/separator/BalancedSeparator.hpp"
#include "techniques/separator/LiptonTarjanFundamentalCycle.hpp"

namespace Koala {

LiptonTarjanPlanarSeparator::LiptonTarjanPlanarSeparator(const NetworKit::Graph &graph,
                                                         const std::vector<double> &costs)
    : BalancedSeparator(graph), vertex_cost(costs) {}

LiptonTarjanPlanarSeparator::LiptonTarjanPlanarSeparator(const NetworKit::Graph &graph)
    : BalancedSeparator(graph), vertex_cost(graph.upperNodeIdBound(), 1.0) {}

void LiptonTarjanPlanarSeparator::run() {
    clean_partitions();

    if (graph.numberOfNodes() == 0) {
        hasRun = true;
        return;
    }

    // normalize total cost to 1
    double totalCost = 0.0;

    assert(totalCost <= 1.0);
    graph.forNodes([&](NetworKit::node v) { totalCost += vertex_cost[v]; });
    if (totalCost > 0.0) {
        graph.forNodes([&](NetworKit::node v) { vertex_cost[v] /= totalCost; });
    }

    // Step 2: Find connected components of the graph G
    auto componentsAlgorithm = NetworKit::ConnectedComponents(graph);
    componentsAlgorithm.run();

    auto components = componentsAlgorithm.getComponents();

    if (are_connected_components_eligible_for_partition(components)) {
        DEBUG("It is possible to find a partition without finding cycle.\n Using "
              "find_separator_from_components");
        find_separator_from_components(components);
    } else {
        // Step 3: Perform BFS
        auto G = componentsAlgorithm.extractLargestConnectedComponent(graph, false);

        // find the root for BFS spanning tree
        NetworKit::node root = NetworKit::none;
        G.forNodes([&](NetworKit::node v) {
            if (root == NetworKit::none)
                root = v;
        });

        auto [lvl, parent] = perform_BFS_and_find_spanning_tree(G, root);
        auto verticesAtLevel = find_number_of_vertices_at_level(lvl);

        auto [l0, l1, l2] = find_partition_levels(G, verticesAtLevel, lvl);

        if (compute_middle_cost(G, lvl, l0, l2) <= 2.0 / 3) {
            DEBUG("Cost of middle part is <= 2/3. Fallback to "
                  "extract_level_partition");
            extract_level_partition(G, lvl, l0, l2);
            hasRun = true;
            return;
        }

        // Step 6: Remove vertices at level >= l2, contract levels <= l0 into
        // vertex x
        auto embedding = PlanarGraphTools::findPlanarEmbedding(graph);
        auto [H, x] = find_contracted_subgraph(G, embedding, lvl, l0, l2);

        // Step 7 compute spanning tree of H and costs of subtrees hanging off
        // of it
        auto [lvlH, parentH] = perform_BFS_and_find_spanning_tree(H, x);

        H = PlanarGraphTools::makeMaximalPlanar(H);
        auto embeddingH = PlanarGraphTools::findPlanarEmbedding(H);
        auto costsH = find_subtree_costs(H, parentH, vertex_cost, x);

        // Steps 8 and 9
        LiptonTarjanFundamentalCycle fundamentalCycleAlgorithm(H, embeddingH, parentH, costsH,
                                                               vertex_cost, x);
        fundamentalCycleAlgorithm.run();
        auto fundamental_cycle_opt = fundamentalCycleAlgorithm.getFundamentalCycle();
        if (!fundamental_cycle_opt) {
            ERROR("Fundamental cycle returned by fundamental cycle algo is "
                  "empty");
            hasRun = true;
            return;
        }

        auto fundamental_cycle = fundamental_cycle_opt.value();

        // ---- Step 10: extract separator and partitions back to G ----
        extract_separator_and_partition_from_cycle(graph, lvl, l0, l2, fundamental_cycle, x);
    }
    hasRun = true;
}

std::tuple<int, int, int>
LiptonTarjanPlanarSeparator::find_partition_levels(const NetworKit::Graph &G,
                                                   std::vector<NetworKit::node> &verticesAtLevel,
                                                   std::vector<int> &lvl) {
    std::vector<NetworKit::count> prefixVertices(verticesAtLevel.size());
    std::partial_sum(verticesAtLevel.begin(), verticesAtLevel.end(), prefixVertices.begin());

    std::vector<double> costAtLevel(verticesAtLevel.size(), 0.0);
    double costG = 0.0;
    G.forNodes([&](NetworKit::node v) {
        costAtLevel[lvl[v]] += vertex_cost[v];
        costG += vertex_cost[v];
    });
    std::vector<double> prefixCost(costAtLevel.size());
    std::partial_sum(costAtLevel.begin(), costAtLevel.end(), prefixCost.begin());

    // Step 4: Find level l1, and k
    // l1 is the smallest level such that the total cost of levels 0..l1 is at
    // least half of the total cost of G.
    NetworKit::node l1 =
        std::min(static_cast<NetworKit::node>(
                     std::lower_bound(prefixCost.begin(), prefixCost.end(), 0.5 * costG)
                     - prefixCost.begin()),
                 static_cast<NetworKit::node>(verticesAtLevel.size() - 1));

    // k is the number of vertices in levels 0..l1 .
    NetworKit::count k = prefixVertices[l1];

    // Step 5: Find levels l0 and l2
    // l0 is the largest level <= l1 such that |L(l0)| + 2*(l1-l0) <= 2*sqrt(k).
    NetworKit::node l0 = 0;
    for (int i = static_cast<int>(l1); i >= 0; i--) {
        double lhs = verticesAtLevel[i] + 2.0 * (static_cast<int>(l1) - i);
        double rhs = 2.0 * std::sqrt(static_cast<double>(k));
        if (lhs <= rhs) {
            l0 = static_cast<NetworKit::node>(i);
            break;
        }
    }
    // l2 is the smallest level > l1 such that
    // |L(l2)| + 2*(l2-l1-1) <= 2*sqrt(n-k)

    int r = static_cast<int>(verticesAtLevel.size()) - 1;
    int l2 = r + 1;
    for (int i = l1 + 1; i <= r + 1; i++) {
        double size_i = (i == r + 1) ? 0.0 : static_cast<double>(verticesAtLevel[i]);
        double lhs = size_i + 2.0 * (i - l1 - 1);
        double rhs = 2.0 * std::sqrt(static_cast<double>(G.numberOfNodes() - k));
        if (lhs <= rhs) {
            l2 = i;
            break;
        }
    }

    return {l0, l1, l2};
}

std::vector<double> LiptonTarjanPlanarSeparator::find_subtree_costs(
    const NetworKit::Graph &H, const std::vector<NetworKit::node> &parentH,
    const std::vector<double> &vertex_cost, NetworKit::node root) {
    std::vector<std::vector<NetworKit::node>> childrenH(H.upperNodeIdBound(),
                                                        std::vector<NetworKit::node>());
    H.forNodes([&](NetworKit::node v) {
        if (parentH[v] != NetworKit::none) {
            childrenH[parentH[v]].push_back(v);
        }
    });
    std::vector<double> costsH(H.upperNodeIdBound(), 0.0);
    std::vector<std::pair<NetworKit::node, size_t>> stk;
    costsH[root] = 0.0;
    stk.emplace_back(root, 0);
    while (!stk.empty()) {
        NetworKit::node v = stk.back().first;
        size_t &idx = stk.back().second;
        if (idx < childrenH[v].size()) {
            NetworKit::node c = childrenH[v][idx++];
            costsH[c] = vertex_cost[c];
            stk.emplace_back(c, 0);
        } else {
            NetworKit::node finished = v;
            stk.pop_back();
            if (!stk.empty())
                costsH[stk.back().first] += costsH[finished];
        }
    }
    return costsH;
}

std::vector<NetworKit::count>
LiptonTarjanPlanarSeparator::find_number_of_vertices_at_level(const std::vector<int> &lvl) {
    int maxLevel = -1;
    for (auto l : lvl) {
        if (l != -1)
            maxLevel = std::max(maxLevel, l);
    }
    std::vector<NetworKit::count> verticesAtLevel(maxLevel + 1, 0);
    for (auto l : lvl) {
        if (l != -1)
            verticesAtLevel[l]++;
    }
    return verticesAtLevel;
}

std::pair<std::vector<int>, std::vector<NetworKit::node>>
LiptonTarjanPlanarSeparator::perform_BFS_and_find_spanning_tree(const NetworKit::Graph &G,
                                                                NetworKit::node startNode) {
    NetworKit::count n = G.upperNodeIdBound();

    std::vector<int> lvl(n, -1);
    std::vector<NetworKit::node> parent(n, NetworKit::none);

    NetworKit::BFS bfs(G, startNode);
    bfs.run();

    G.forNodes([&](NetworKit::node v) {
        if (bfs.getDistances()[v] < std::numeric_limits<double>::max()) {
            lvl[v] = static_cast<NetworKit::node>(bfs.getDistances()[v]);
            auto preds = bfs.getPredecessors(v);
            if (!preds.empty()) {
                parent[v] = preds[0];
            }
        }
    });

    return {lvl, parent};
}

std::pair<NetworKit::Graph, NetworKit::node>
LiptonTarjanPlanarSeparator::find_contracted_subgraph(const NetworKit::Graph &G,
                                                      planar_embedding_t &embedding,
                                                      std::vector<int> &lvl, int l0, int l2) {
    std::unordered_set<NetworKit::node> verticesToKeep;
    G.forNodes([&](NetworKit::node v) {
        if (lvl[v] < l2)
            verticesToKeep.insert(v);
    });

    auto H = NetworKit::GraphTools::subgraphFromNodes(G, verticesToKeep);
    NetworKit::node x = H.addNode();

    std::vector<bool> isConnectedToX(H.upperNodeIdBound(), false);

    // Walk around the subtree (levels 0..l0) using the embedding
    // and rewire boundary edges to x
    std::unordered_set<NetworKit::node> subtreeNodes;
    G.forNodes([&](NetworKit::node v) {
        if (lvl[v] <= l0)
            subtreeNodes.insert(v);
    });

    for (auto v : subtreeNodes) {
        for (auto w : embedding[v]) {
            if (subtreeNodes.count(w) > 0)
                continue; // both in subtree, skip
            if (!H.hasNode(w))
                continue; // deleted (level >= l2), skip
            if (!isConnectedToX[w]) {
                isConnectedToX[w] = true;
                H.addEdge(x, w);
            }
        }
    }

    // Remove the contracted subtree nodes from H
    for (auto v : subtreeNodes) {
        H.removeNode(v);
    }

    return {H, x};
}

bool LiptonTarjanPlanarSeparator::are_connected_components_eligible_for_partition(
    std::vector<std::vector<NetworKit::node>> &components) {
    double highestCost = 0.0;
    for (const auto &cc : components) {
        double ccCost = 0.0;
        for (auto v : cc)
            ccCost += vertex_cost[v];
        highestCost = std::max(highestCost, ccCost);
    }

    return highestCost <= 2.0 / 3.0;
}

void LiptonTarjanPlanarSeparator::find_separator_from_components(
    std::vector<std::vector<NetworKit::node>> &components) {
    double mostExpensiveComponentCost = 0.0;
    size_t mostExpensiveComponentId = 0;

    auto findComponentCost = [&](std::vector<NetworKit::node> &component) -> double {
        double cost = 0.0;
        for (auto node : component) {
            cost += vertex_cost[node];
        }

        return cost;
    };

    // find most expensive component
    for (size_t i = 0; i < components.size(); i++) {
        double cost = findComponentCost(components[i]);
        if (cost > mostExpensiveComponentCost) {
            mostExpensiveComponentCost = cost;
            mostExpensiveComponentId = i;
        }
    }

    if (mostExpensiveComponentCost > 1.0 / 3 && mostExpensiveComponentCost <= 2.0 / 3) {
        partition.A.insert(partition.A.end(), components[mostExpensiveComponentId].begin(),
                           components[mostExpensiveComponentId].end());

        for (size_t i = 0; i < components.size(); i++) {
            if (i != mostExpensiveComponentId) {
                partition.B.insert(partition.B.end(), components[i].begin(), components[i].end());
            }
        }
    } else {
        size_t thresholdId = 0;
        double totalCost = 0.0;

        for (; thresholdId < components.size(); thresholdId++) {
            totalCost += findComponentCost(components[thresholdId]);
            if (totalCost > 1.0 / 3)
                break;
        }

        for (size_t i = 0; i < components.size(); i++) {
            if (i <= thresholdId)
                partition.A.insert(partition.A.end(), components[i].begin(), components[i].end());
            else
                partition.B.insert(partition.B.end(), components[i].begin(), components[i].end());
        }
    }
}

void LiptonTarjanPlanarSeparator::extract_separator_and_partition_from_cycle(
    const NetworKit::Graph &G, const std::vector<int> &lvl, int l0, int l2,
    std::vector<NetworKit::node> &finalCycle, NetworKit::node x) {
    NetworKit::Graph graphCopy = G;
    std::unordered_set<NetworKit::node> uniqueSeparator;

    graphCopy.forNodes([&](NetworKit::node v) {
        if (lvl[v] == l0 || lvl[v] == l2)
            uniqueSeparator.insert(v);
    });

    for (auto v : finalCycle) {
        if (v != x)
            uniqueSeparator.insert(v);
    }

    partition.separator.assign(uniqueSeparator.begin(), uniqueSeparator.end());

    for (auto v : partition.separator) {
        graphCopy.removeNode(v);
    }

    NetworKit::ConnectedComponents connectedComponentsAlgo =
        NetworKit::ConnectedComponents(graphCopy);
    connectedComponentsAlgo.run();
    auto connectedComponents = connectedComponentsAlgo.getComponents();
    find_separator_from_components(connectedComponents);
}

void LiptonTarjanPlanarSeparator::extract_level_partition(const NetworKit::Graph &G,
                                                          const std::vector<int> &lvl, int l0,
                                                          int l2) {
    struct Part {
        int id;
        double cost;
        std::vector<NetworKit::node> nodes;
    };

    int maxLvl = 0;
    for (auto l : lvl) {
        if (l != -1)
            maxLvl = std::max(maxLvl, l);
    }

    auto findPartCost = [&](int id, int lower, int upper, bool empty = false) -> Part {
        double cost = 0.0;
        std::vector<NetworKit::node> partVertices;

        if (!empty) {
            G.forNodes([&](NetworKit::node t) {
                if (lvl[t] != -1 && lvl[t] >= lower && lvl[t] <= upper) {
                    cost += vertex_cost[t];
                    partVertices.push_back(t);
                }
            });
        }
        return {id, cost, partVertices};
    };

    bool part1Empty = (l0 == 0);
    int upper0 = part1Empty ? 0 : l0 - 1;

    std::vector<Part> parts({findPartCost(1, 0, upper0, part1Empty),
                             findPartCost(2, l0 + 1, l2 - 1, (l0 + 1 > l2 - 1)),
                             findPartCost(3, l2 + 1, maxLvl, (l2 + 1 > maxLvl))});
    auto mostExpensivePart = *(std::ranges::max_element(parts, {}, &Part::cost));

    partition.A.insert(partition.A.end(), mostExpensivePart.nodes.begin(),
                       mostExpensivePart.nodes.end());
    for (auto &part : parts) {
        if (part.id != mostExpensivePart.id) {
            partition.B.insert(partition.B.end(), part.nodes.begin(), part.nodes.end());
        }
    }
    G.forNodes([&](NetworKit::node t) {
        if (lvl[t] == l0 || lvl[t] == l2) {
            partition.separator.push_back(t);
        }
    });
}

double LiptonTarjanPlanarSeparator::compute_middle_cost(const NetworKit::Graph &G,
                                                        const std::vector<int> &lvl, int l0,
                                                        int l2) {
    double middleCost = 0.0;
    G.forNodes([&](NetworKit::node t) {
        if (lvl[t] >= l0 + 1 && lvl[t] <= l2 - 1)
            middleCost += vertex_cost[t];
    });

    return middleCost;
}

void LiptonTarjanPlanarSeparator::clean_partitions() {
    partition.separator.clear();
    partition.A.clear();
    partition.B.clear();
}

} // namespace Koala
