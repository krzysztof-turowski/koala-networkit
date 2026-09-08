#include "techniques/separator/LiptonTarjanFundamentalCycle.hpp"

#include <algorithm>
#include <optional>
#include <unordered_set>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/EdgeUtils.hpp>

#include "shortest_path/planar/SuitableRDivision.hpp"

namespace {
bool isTreeEdge(NetworKit::node u, NetworKit::node v,
                std::vector<NetworKit::node>& parent) {
    return u == parent[v] || v == parent[u];
}

std::pair<NetworKit::node, NetworKit::node> findNonTreeEdge(
    NetworKit::Graph& H, std::vector<NetworKit::node>& parent) {
    NetworKit::node v1 = NetworKit::none;
    NetworKit::node w1 = NetworKit::none;

    // find non tree edge (v1, w1)
    H.forEdges([&](NetworKit::node v, NetworKit::node w) {
        if (v1 != NetworKit::none) return;  // already found
        if (parent[v] != w && parent[w] != v) {
            v1 = v;
            w1 = w;
        }
    });

    return {v1, w1};
}

bool isOnInsideArc(int pos, int posNext, int posPrev) {
    if (posNext < posPrev)
        return pos > posNext && pos < posPrev;
    else
        return pos > posNext || pos < posPrev;
}
}  // namespace

namespace Koala {
LiptonTarjanFundamentalCycle::LiptonTarjanFundamentalCycle(
    NetworKit::Graph& G, planar_embedding_t& embedding,
    std::vector<NetworKit::node>& parent, std::vector<double>& costsOfSubtree,
    std::vector<double>& vertexCost, NetworKit::node root)
    : G(G),
      embedding(embedding),
      parent(parent),
      costs_of_subtree(costsOfSubtree),
      vertex_cost(vertexCost),
      root(root) {}

std::optional<cycle_t> LiptonTarjanFundamentalCycle::getFundamentalCycle() {
    assureFinished();
    return fundamental_cycle;
}

void LiptonTarjanFundamentalCycle::run() {
    auto [v1, w1] = findNonTreeEdge(G, parent);
    find_indexing_map();

    auto initial_cycle_opt = build_fundamental_cycle(v1, w1, parent);

    if (!initial_cycle_opt) {
        ERROR("Fundamental cycle not found by build_fundamental_cycle");
        fundamental_cycle = std::nullopt;
        hasRun = true;
        return;
    }

    auto initial_cycle = initial_cycle_opt.value();
    if (initial_cycle.size() < 3) {
        ERROR(
            "Fundamental cycle found by build_fundamental_cycle is degenerate");
        fundamental_cycle = std::nullopt;
        hasRun = true;
        return;
    }

    auto sides_cost = compute_sides_initial_cost(initial_cycle);

    auto final_cycle = shrink_fundamental_cycle(initial_cycle, sides_cost,
                                                idx_of, root, v1, w1);
    if (final_cycle.empty()) {
        ERROR("Shrink fundamental cycle returned an empty cycle");
        fundamental_cycle = std::nullopt;
        hasRun = true;
        return;
    }

    fundamental_cycle = std::move(final_cycle);

    hasRun = true;
}

std::optional<cycle_t> LiptonTarjanFundamentalCycle::build_fundamental_cycle(
    NetworKit::node v1, NetworKit::node w1,
    std::vector<NetworKit::node> parent) {
    // Degenerate H (no non-tree edge): fall back to a level-based separator.
    if (v1 == NetworKit::none) {
        return std::nullopt;
    }

    // Collect ancestors of v1 and w1
    std::vector<NetworKit::node> pathV, pathW;

    // Walk v1 to root
    for (auto v = v1; v != NetworKit::none; v = parent[v]) pathV.push_back(v);

    // Walk w1 to root
    for (auto w = w1; w != NetworKit::none; w = parent[w]) pathW.push_back(w);

    // Find LCA — first common vertex
    std::unordered_set<NetworKit::node> ancestorsV(pathV.begin(), pathV.end());
    NetworKit::node lca = NetworKit::none;
    for (auto w : pathW) {
        if (ancestorsV.count(w)) {
            lca = w;
            break;
        }
    }

    // build cycle
    cycle_t cycle;
    for (auto v : pathV) {
        cycle.push_back(v);
        if (v == lca) break;
    }

    std::vector<NetworKit::node> pathWToLCA;
    for (auto w : pathW) {
        if (w == lca) break;
        pathWToLCA.push_back(w);
    }
    std::reverse(pathWToLCA.begin(), pathWToLCA.end());
    cycle.insert(cycle.end(), pathWToLCA.begin(), pathWToLCA.end());

    return cycle;
}

CostComputation LiptonTarjanFundamentalCycle::compute_sides_initial_cost(
    const cycle_t& cycle) {
    double arcTrueCost = 0.0;
    double arcFalseCost = 0.0;

    int cycleSize = cycle.size();

    for (int i = 0; i < cycleSize; i++) {
        auto prev = find_node_in_cycle_at_ith_position(cycle, i - 1);
        auto v = find_node_in_cycle_at_ith_position(cycle, i);
        auto u = find_node_in_cycle_at_ith_position(cycle, i + 1);

        auto neighborsV = embedding[v];
        int posU = idx_of[v][u];
        int posPrev = idx_of[v][prev];

        for (int j = 0; j < static_cast<int>(neighborsV.size()); j++) {
            auto w = neighborsV[j];
            if (w == prev || w == u) continue;

            double cost = 0.0;
            if (parent[w] == v) {
                cost = costs_of_subtree[w];
            } else if (parent[v] == w) {
                cost = costs_of_subtree[root] - costs_of_subtree[v];
            } else {
                continue;
            }

            if (isOnInsideArc(j, posU, posPrev)) {
                arcTrueCost += cost;
            } else {
                arcFalseCost += cost;
            }
        }
    }

    double insideCost = fmax(arcTrueCost, arcFalseCost);
    double outsideCost = fmin(arcTrueCost, arcFalseCost);
    bool isTrueArcInside = arcTrueCost >= arcFalseCost;

    return {insideCost, outsideCost, isTrueArcInside};
}

void LiptonTarjanFundamentalCycle::find_indexing_map() {
    std::vector<std::unordered_map<NetworKit::node, int>> idx_map(
        G.upperNodeIdBound());

    G.forNodes([&](NetworKit::node u) {
        const auto& rotation = embedding[u];

        for (int i = 0; i < static_cast<int>(rotation.size()); ++i) {
            idx_map[u][rotation[i]] = i;
        }
    });

    idx_of = idx_map;
}

NetworKit::node LiptonTarjanFundamentalCycle::getApex(NetworKit::node v,
                                                      NetworKit::node u,
                                                      int dir) {
    const auto& neighborsV = embedding[v];
    int deg = static_cast<int>(neighborsV.size());

    int j = idx_of[v].at(u);
    return neighborsV[wrapIndex(j + dir, deg)];
}

std::tuple<std::vector<NetworKit::node>, std::vector<bool>, NetworKit::node>
LiptonTarjanFundamentalCycle::compute_path_to_cycle(
    NetworKit::node apex, NetworKit::node curVi, NetworKit::node curWi,
    std::vector<bool>& isOnCycle) {
    std::vector<NetworKit::node> Ppath;
    std::vector<bool> isOnPath(G.upperNodeIdBound(), false);
    NetworKit::node z = apex;
    while (z != NetworKit::none && isOnCycle[z] == false) {
        Ppath.push_back(z);
        isOnPath[z] = true;
        z = parent[z];
    }

    // it is possible that the root is inside of boundary, as we deduce inside
    // vs outside based on cost, not on topology. For that we need to find LCA
    // of curWi and curVi on the cycle and compute the path as concat of y ~>
    // LCA(y, z) and z ~> LCA(y,z)
    if (z == NetworKit::none) {
        auto getLCA = [&](NetworKit::node a, NetworKit::node b) {
            std::unordered_set<NetworKit::node> ancA;
            for (auto c = a; c != NetworKit::none; c = parent[c])
                ancA.insert(c);
            for (auto c = b; c != NetworKit::none; c = parent[c]) {
                if (ancA.count(c)) return c;
            }
            return NetworKit::none;
        };

        z = getLCA(curVi, curWi);

        std::vector<NetworKit::node> pY;
        for (auto c = apex; c != NetworKit::none; c = parent[c])
            pY.push_back(c);

        std::vector<NetworKit::node> pZ;
        for (auto c = z; c != NetworKit::none; c = parent[c]) pZ.push_back(c);

        std::unordered_set<NetworKit::node> setY(pY.begin(), pY.end());
        NetworKit::node aCommon = NetworKit::none;
        for (auto c : pZ) {
            if (setY.count(c)) {
                aCommon = c;
                break;
            }
        }
        aCommon = root;

        std::vector<NetworKit::node> fullSeq;
        for (auto c : pY) {
            fullSeq.push_back(c);
            if (c == aCommon) break;
        }
        std::vector<NetworKit::node> pathDown;
        for (auto c : pZ) {
            if (c == aCommon) break;
            pathDown.push_back(c);
        }
        std::reverse(pathDown.begin(), pathDown.end());
        for (auto c : pathDown) fullSeq.push_back(c);

        for (auto c : Ppath) isOnPath[c] = false;
        Ppath.clear();
        if (!fullSeq.empty()) {
            for (size_t i = 0; i < fullSeq.size() - 1; i++) {
                Ppath.push_back(fullSeq[i]);
                isOnPath[fullSeq[i]] = true;
            }
        }
    }

    return {Ppath, isOnPath, z};
}

NetworKit::node
LiptonTarjanFundamentalCycle::find_node_in_cycle_at_ith_position(
    const cycle_t& cycle, int i) {
    return cycle[wrapIndex(i, static_cast<int>(cycle.size()))];
}
int LiptonTarjanFundamentalCycle::wrapIndex(int i, int size) {
    return ((i % size) + size) % size;
}

cycle_t LiptonTarjanFundamentalCycle::shrink_fundamental_cycle(
    const cycle_t& initial_cycle, const CostComputation& sides_cost,
    const index_map_t& idx_of, NetworKit::node root, NetworKit::node v1,
    NetworKit::node w1) {
    std::vector<bool> isOnCycle(G.upperNodeIdBound(), false);
    for (auto v : initial_cycle) {
        isOnCycle[v] = true;
    }

    std::vector<NetworKit::node> cyclePrev(G.upperNodeIdBound(),
                                           NetworKit::none);
    std::vector<NetworKit::node> cycleNext(G.upperNodeIdBound(),
                                           NetworKit::none);

    for (size_t i = 0; i < initial_cycle.size(); i++) {
        NetworKit::node cur = initial_cycle[i];
        NetworKit::node prev =
            initial_cycle[wrapIndex(i - 1, initial_cycle.size())];
        NetworKit::node next =
            initial_cycle[wrapIndex(i + 1, initial_cycle.size())];
        cycleNext[cur] = next;
        cyclePrev[cur] = prev;
    }

    int dir = sides_cost.insideIsClockwiseArc ? -1 : 1;

    NetworKit::count nodeBound = G.upperNodeIdBound();
    double threshold = 2.0 / 3.0 * costs_of_subtree[root];

    NetworKit::node curVi = v1;
    NetworKit::node curWi = w1;
    double insideCost = sides_cost.insideCost;

    // safeguard iterGuard
    long iterGuard = 0;
    const long iterCap = 4 * static_cast<long>(nodeBound);

    while (insideCost > threshold) {
        if (++iterGuard > iterCap) {
            ERROR("Iteration cap reached");
            return {};
        }

        NetworKit::node y = getApex(curVi, curWi, dir);

        if (isTreeEdge(curVi, y, parent) || isTreeEdge(curWi, y, parent)) {
            auto insertBetween = [&](NetworKit::node u, NetworKit::node v,
                                     NetworKit::node mid) {
                if (cycleNext[u] == v) {
                    cycleNext[u] = mid;
                    cyclePrev[mid] = u;
                    cycleNext[mid] = v;
                    cyclePrev[v] = mid;
                } else {
                    cyclePrev[u] = mid;
                    cycleNext[mid] = u;
                    cyclePrev[mid] = v;
                    cycleNext[v] = mid;
                }
            };

            auto bypassNode = [&](NetworKit::node u, NetworKit::node old,
                                  NetworKit::node v) {
                if (cycleNext[u] == old)
                    cycleNext[u] = v;
                else
                    cyclePrev[u] = v;
                if (cycleNext[v] == old)
                    cycleNext[v] = u;
                else
                    cyclePrev[v] = u;
            };

            bool edgeViY_is_tree = isTreeEdge(curVi, y, parent);

            if (!isOnCycle[y]) {
                insideCost -= vertex_cost[y];
                isOnCycle[y] = true;

                insertBetween(curVi, curWi, y);

                if (edgeViY_is_tree)
                    curVi = y;
                else
                    curWi = y;

            } else {
                if (edgeViY_is_tree) {
                    isOnCycle[curVi] = false;
                    bypassNode(y, curVi, curWi);
                    curVi = y;
                } else {
                    isOnCycle[curWi] = false;
                    bypassNode(y, curWi, curVi);
                    curWi = y;
                }
            }
        } else {
            auto [Ppath, isOnPath, z] =
                compute_path_to_cycle(y, curVi, curWi, isOnCycle);
            double PpathCost = 0.0;
            for (auto v : Ppath) {
                PpathCost += vertex_cost[v];
            }

            struct EdgeScanner {
                NetworKit::node beg, cur = beg, prev;
                int embCur = -1, embStop = -1;
                int rotDir;
                double cost = 0.0;
                bool forward = false, done = false, isInitRequired = true;

                EdgeScanner(NetworKit::node beg, NetworKit::node prev,
                            int rotDir, std::vector<NetworKit::node>& cyclePrev)
                    : beg(beg), prev(prev), rotDir(rotDir) {
                    forward = cyclePrev[cur] == prev;
                }
            };

            auto nextInCycleDir = [&](NetworKit::node u,
                                      bool forward) -> NetworKit::node {
                return forward ? cycleNext[u] : cyclePrev[u];
            };

            EdgeScanner scannerV{curVi, curWi, dir, cyclePrev};
            EdgeScanner scannerW{curWi, curVi, -dir, cyclePrev};

            auto walkThePath = [&](EdgeScanner& scanner) {
                for (int i = static_cast<int>(Ppath.size()) - 1; i >= 0; i--) {
                    NetworKit::node pathPrev =
                        (i == static_cast<int>(Ppath.size() - 1))
                            ? z
                            : Ppath[i + 1];
                    NetworKit::node pathNext =
                        (i == 0) ? scanner.beg : Ppath[i - 1];
                    NetworKit::node pathCur = Ppath[i];
                    int deg = static_cast<int>(embedding[pathCur].size());

                    int embBeg = wrapIndex(
                        idx_of[pathCur].at(pathPrev) + scanner.rotDir, deg);
                    int embEnd = idx_of[pathCur].at(pathNext);

                    for (int j = embBeg; j != embEnd;
                         j = wrapIndex(j + scanner.rotDir, deg)) {
                        NetworKit::node w = embedding[pathCur][j];
                        if (parent[w] == pathCur) {
                            scanner.cost += costs_of_subtree[w];
                        } else if (parent[pathCur] == w) {
                            scanner.cost += costs_of_subtree[root] -
                                            costs_of_subtree[pathCur];
                        }
                    }
                }
            };

            // a one step at at time boundary walking strucutre to ensure O(n)
            // time as described in Step 9 of the algorithm found from Lipton,
            // Richard J., and Robert Endre Tarjan. "A separator theorem for
            // planar graphs." SIAM Journal on Applied Mathematics 36.2 (1979):
            // 177-189.
            auto stepBoundary = [&](EdgeScanner& scanner) {
                if (scanner.isInitRequired) {
                    fflush(stdout);

                    scanner.embCur = wrapIndex(
                        idx_of[scanner.cur].at(scanner.prev) + scanner.rotDir,
                        embedding[scanner.cur].size());

                    NetworKit::node back =
                        Ppath.empty() ? scanner.beg : Ppath.back();
                    NetworKit::node next =
                        scanner.cur == z
                            ? back
                            : nextInCycleDir(scanner.cur, scanner.forward);
                    scanner.embStop = wrapIndex(idx_of[scanner.cur].at(next),
                                                embedding[scanner.cur].size());
                    scanner.isInitRequired = false;
                }

                if (scanner.embCur == scanner.embStop) {
                    if (scanner.cur == z) {
                        scanner.done = true;
                        walkThePath(scanner);
                        return;
                    }

                    scanner.prev = scanner.cur;
                    scanner.cur = nextInCycleDir(scanner.cur, scanner.forward);
                    scanner.isInitRequired = true;
                } else {
                    NetworKit::node w = embedding[scanner.cur][scanner.embCur];
                    if (parent[w] == scanner.cur) {
                        scanner.cost += costs_of_subtree[w];
                    } else if (parent[scanner.cur] == w) {
                        scanner.cost += costs_of_subtree[root] -
                                        costs_of_subtree[scanner.cur];
                    }

                    scanner.embCur = wrapIndex(scanner.embCur + scanner.rotDir,
                                               embedding[scanner.cur].size());
                }
            };

            while (!scannerV.done && !scannerW.done) {
                stepBoundary(scannerV);
                stepBoundary(scannerW);
            }

            EdgeScanner finishedScanner = scannerV.done ? scannerV : scannerW;
            EdgeScanner unfinishedScanner = scannerV.done ? scannerW : scannerV;
            const double finishedSideCost = finishedScanner.cost;
            const double unfinishedSideCost =
                insideCost - finishedSideCost - PpathCost;
            const bool keepFinishedSide =
                finishedSideCost >= unfinishedSideCost;
            EdgeScanner keptScanner =
                keepFinishedSide ? finishedScanner : unfinishedScanner;

            insideCost = std::max(finishedSideCost, unfinishedSideCost);

            for (auto v : Ppath) {
                isOnCycle[v] = true;
            }

            // rewire cycle lists
            auto& dirPrev = keptScanner.forward ? cyclePrev : cycleNext;
            auto& dirNext = keptScanner.forward ? cycleNext : cyclePrev;

            dirPrev[keptScanner.beg] = y;
            dirNext[z] = Ppath.empty() ? keptScanner.beg : Ppath.back();

            for (int i = static_cast<int>(Ppath.size()) - 1; i >= 0; i--) {
                NetworKit::node pathPrev =
                    (i == static_cast<int>(Ppath.size() - 1)) ? z
                                                              : Ppath[i + 1];
                NetworKit::node pathNext =
                    (i == 0) ? keptScanner.beg : Ppath[i - 1];
                NetworKit::node pathCur = Ppath[i];

                dirPrev[pathCur] = pathPrev;
                dirNext[pathCur] = pathNext;
            }

            if (keptScanner.beg == curVi) {
                curWi = y;
            } else {
                curVi = y;
            }
        }
    }
    cycle_t finalCycle;
    NetworKit::node curr = curVi;
    do {
        finalCycle.push_back(curr);
        curr = cycleNext[curr];
    } while (curr != curVi);

    return finalCycle;
}

}  // namespace Koala
