#include "planar_sssp/PlanarSSSP.hpp"

namespace Koala {

PlanarSSSP::PlanarSSSP(NetworKit::Graph& graph, NetworKit::node source, NetworKit::node target)
    : graph(graph), source(source), target(target), distances(graph.numberOfNodes()) {}

NetworKit::edgeweight PlanarSSSP::distance(NetworKit::node t) {
    assureFinished();
    return distances[t];
}

const std::vector<NetworKit::edgeweight> PlanarSSSP::getDistances() {
    assureFinished();
    return distances;
}

} /* namespace Koala */
