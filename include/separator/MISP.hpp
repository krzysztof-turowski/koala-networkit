#include "networkit/Globals.hpp"
#include "networkit/base/Algorithm.hpp"
#include "networkit/graph/Graph.hpp"
#include <functional>
#include <map>
#include <vector>
namespace Koala {
template <typename T> class MISP : public NetworKit::Algorithm {
public:
  using ComponentSolver =
      std::function<std::vector<T>(const NetworKit::Graph &component)>;
  void run() override;

  explicit MISP(
      const NetworKit::Graph &G, double epsilon,
      ComponentSolver componentSolver,
      std::optional<std::map<NetworKit::node, double>> costs = std::nullopt);
  std::vector<T> maximum_independent_set;

private:
  NetworKit::Graph graph;
  double epsilon;
  ComponentSolver componentSolver;
  std::optional<std::map<NetworKit::node, double>> costs;
};

} // namespace Koala
