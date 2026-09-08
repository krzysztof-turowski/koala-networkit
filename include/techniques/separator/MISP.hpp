#pragma once
#include <functional>
#include <queue>
#include <vector>
#include "networkit/Globals.hpp"
#include "networkit/base/Algorithm.hpp"
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {
template <typename T>
class MISP : public NetworKit::Algorithm {
public:
    using ComponentSolver = std::function<std::vector<T>(const NetworKit::Graph &component)>;
    void run() override;

    explicit MISP(const NetworKit::Graph &G, double epsilon, ComponentSolver componentSolver,
                  std::optional<std::vector<double>> costs = std::nullopt);

    std::vector<T> getIndependentSet();

private:
    std::vector<T> maximum_independent_set;

    std::vector<NetworKit::node> find_espilon_planar_separator();
    double get_component_cost(const std::vector<NetworKit::node> &connected_component);
    void
    processSide(const std::vector<NetworKit::node> &subGraph,
                std::unordered_map<NetworKit::node, std::vector<NetworKit::node>> &idToComponentMap,
                std::queue<int> &Q, int &idCnt);

    NetworKit::Graph graph;
    double epsilon;
    ComponentSolver componentSolver;
    std::vector<double> vertex_cost;
};
template class MISP<NetworKit::Edge>;
template class MISP<NetworKit::node>;

} // namespace Koala
