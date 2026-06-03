#include "flow/GoldbergTarjanPushRelabelMaximumFlow.hpp"

namespace Koala {

NetworKit::node GoldbergTarjanPushRelabelMaximumFlow::get_active_vertex() {
    for (NetworKit::node v = 0; v < graph.upperNodeIdBound(); ++v) {
        if (v != source && v != target && graph.hasNode(v) && excess[v] > 0) {
            return v;
        }
    }
    return NetworKit::none;
}

NetworKit::node GoldbergTarjanPushRelabelMaximumFlow::get_admissible_residual_edge(
        NetworKit::node v) {
    NetworKit::node result = NetworKit::none;
    graph.forNeighborsOf(v, [&](NetworKit::node u) {
        const NetworKit::Edge e(v, u);
        if (result == NetworKit::none && capacity[e] - flow[e] > 0
                && distance[v] == distance[u] + 1) {
            result = u;
        }
    });
    return result;
}

}  // namespace Koala
