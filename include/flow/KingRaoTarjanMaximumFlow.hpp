#include <flow/MaximumFlow.hpp>

namespace Koala {

/**
 * @ingroup flow
 * The class for the King-Rao-Tarjan maximum flow algorithm
 */
class KingRaoTarjanMaximumFlow final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;

    /**
     * Execute the King-Rao-Tarjan maximum flow algorithm.
     */
    void run();
    int get_flow(const std::pair<NetworKit::node, NetworKit::node>&);

 private:
    std::unordered_map<std::pair<NetworKit::node, NetworKit::node>, int, pair_hash> flow;
    std::map<std::pair<NetworKit::node, NetworKit::node>, int> capacity;
    std::map<NetworKit::node, int> d, excess, hidden_excess;
    std::set<int> positive_excess;
    std::set<std::pair<NetworKit::node, NetworKit::node>> E_star;

    DynamicTree dynamic_tree;
    KRTEdgeDesignator edge_designator;

    int get_visible_excess(NetworKit::node);
    NetworKit::node get_positive_excess_node();
    void update_positive_excess(NetworKit::node);

    void set_flow(const std::pair<NetworKit::node, NetworKit::node>&, int);
    void saturate(const std::pair<NetworKit::node, NetworKit::node>&);
    void add_edge(const std::pair<NetworKit::node, NetworKit::node>&);
    void cut(const std::pair<NetworKit::node, NetworKit::node>&);

    void initialize();
    std::vector<std::pair<NetworKit::node, NetworKit::node>> get_edges_list();
    void tree_push(NetworKit::node, NetworKit::node);
    void relabel(NetworKit::node);
};

}  /* namespace Koala */
