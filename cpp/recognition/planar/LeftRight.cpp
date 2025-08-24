#include <iostream>
#include <list>
#include <set>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>
#include <recognition/planar/PlanarGraphRecognition.hpp>

namespace Koala {
    /// LeftRightPlanarity impl

    bool LeftRightPlanarity::isPlanar(const Graph &G) {
        if (G.numberOfEdges() < 9 || G.numberOfNodes() < 5) return true;
        if (G.numberOfEdges() > 3 * G.numberOfNodes() - 6) return false;

        for (int v = 0; v <= G.numberOfEdges() / G.numberOfNodes(); v++) {
            std::vector<int> dfsHeight2(G.numberOfNodes(), -1);
            dfsHeight2[v] = 0;
            if (!lrAlgorithm(G, v, dfsHeight2)) return false;
        }

        return true;
    }

    bool LeftRightPlanarity::lrAlgorithm(const Graph &G, node root, std::vector<int> &dfsHeight) {
        std::vector<std::vector<LeftRightPlanarity::Fringe>> fringes = {{}};

        std::vector<std::vector<node>> neighbors(G.numberOfNodes());
        for (auto v: G.nodeRange()) {
            for (auto u: G.neighborRange(v)) {
                neighbors[v].push_back(u);
            }
        }

        std::vector<std::pair<node, int>> dfsStack;

        dfsStack.emplace_back(root, 0);

        while (!dfsStack.empty()) {
            auto&[x, pos] = dfsStack.back();
            if (pos < neighbors[x].size()) {
                node y = neighbors[x][pos];

                dfsStack.back().second++;

                if (dfsHeight[y] < 0) {
                    fringes.emplace_back();
                    dfsHeight[y] = dfsHeight[x] + 1;
                    dfsStack.emplace_back(y, 0);
                } else if (dfsHeight[y] < dfsHeight[x]) {
                    fringes.back().emplace_back(dfsHeight[y]);
                }
            } else {
                dfsStack.pop_back();
                if (fringes.size() > 1) {
                    try {
                        mergeFringes(fringes, dfsHeight[dfsStack.back().first]);
                    } catch (...) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    void
    LeftRightPlanarity::mergeFringes(std::vector<std::vector<LeftRightPlanarity::Fringe>> &fringes, int dfs_height) {
        LeftRightPlanarity::Fringe mf = getMergedFringe(fringes.back());
        fringes.pop_back();
        if (!mf.fops.empty()) {
            mf.prune(dfs_height);
            if (!mf.fops.empty()) {
                fringes.back().push_back(mf);
            }
        }
    }

    LeftRightPlanarity::Fringe
    LeftRightPlanarity::getMergedFringe(std::vector<LeftRightPlanarity::Fringe> &upperFringes) {
        if (upperFringes.empty()) return {};
        std::sort(upperFringes.begin(), upperFringes.end());
        LeftRightPlanarity::Fringe result = std::move(upperFringes[0]);
        for (size_t i = 1; i < upperFringes.size(); ++i) {
            result.merge(upperFringes[i]);
        }
        return result;
    }

    /// Fringe impl

    LeftRightPlanarity::FringeOpposedSubset &LeftRightPlanarity::Fringe::H() { return fops.front(); }

    LeftRightPlanarity::FringeOpposedSubset &LeftRightPlanarity::Fringe::L() { return fops.back(); }

    bool LeftRightPlanarity::Fringe::operator<(const LeftRightPlanarity::Fringe &other) const {
        int diff = L().l_lo() - other.L().l_lo();
        return diff != 0 ? diff < 0 : H().l_hi() < other.H().l_hi();
    }

    void LeftRightPlanarity::Fringe::merge(LeftRightPlanarity::Fringe &other) {
        other.mergeTAlikeEdges();
        mergeTOppositeEdgesInto(other);
        if (!H().hasRight()) {
            other.alignDuplicates(L().l_hi());
        } else {
            makeOnionStructure(other);
        }
        if (other.H().hasLeft()) {
            fops.push_front(other.H());
        }
    }

    void LeftRightPlanarity::Fringe::prune(int dfs_height) {
        auto[left_, right_] = lrCondition(dfs_height);
        while (!fops.empty() && (left_ || right_)) {
            if (left_) H().left.pop_front();
            if (right_) H().right.pop_front();
            if (H().left.empty() && H().right.empty()) {
                fops.pop_front();
            } else {
                swapSide();
            }
            if (!fops.empty())
                std::tie(left_, right_) = lrCondition(dfs_height);
        }
    }

    std::pair<bool, bool> LeftRightPlanarity::Fringe::lrCondition(int dfs_height) {
        return {
                H().hasLeft() && H().l_hi() >= dfs_height,
                H().hasRight() && H().r_hi() >= dfs_height
        };
    }

    void LeftRightPlanarity::Fringe::swapSide() {
        if (H().left.empty() || (!H().right.empty() && H().l_lo() > H().r_lo())) {
            std::swap(H().left, H().right);
        }
    }

    void LeftRightPlanarity::Fringe::mergeTAlikeEdges() {
        if (H().hasRight()) {
            throw std::runtime_error("Non-planar");
        }
        for (size_t i = 1; i < fops.size(); ++i) {
            if (!fops[i].right.empty()) throw std::runtime_error("Non-planar");
            H().left.insert(H().left.end(), fops[i].left.begin(), fops[i].left.end());
        }
        fops.erase(fops.begin() + 1, fops.end());
    }

    void LeftRightPlanarity::Fringe::mergeTOppositeEdgesInto(LeftRightPlanarity::Fringe &other) {
        while (!H().hasRight() && !fops.empty() && H().l_hi() > other.H().l_lo()) {
            other.H().right.insert(other.H().right.end(), H().left.begin(), H().left.end());
            fops.pop_front();
        }
    }

    void LeftRightPlanarity::Fringe::alignDuplicates(int dfs_h) {
        if (H().l_lo() == dfs_h) {
            H().left.pop_back();
            swapSide();
        }
    }

    void LeftRightPlanarity::Fringe::makeOnionStructure(LeftRightPlanarity::Fringe &other) {
        int lo = (H().l_hi() < H().r_hi()) ? 0 : 1;
        int hi = 1 - lo;
        auto &lowDeque = lo == 0 ? H().left : H().right;
        auto &highDeque = hi == 0 ? H().left : H().right;
        if (other.H().l_lo() < lowDeque.front()) throw std::runtime_error("Non-planar");
        else if (other.H().l_lo() < highDeque.front()) {
            lowDeque.insert(lowDeque.begin(), other.H().left.rbegin(), other.H().left.rend());
            highDeque.insert(highDeque.begin(), other.H().right.rbegin(), other.H().right.rend());
            other.H().left.clear();
            other.H().right.clear();
        }
    }
}