#include "planar_sssp/PlanarSSSP.hpp"

namespace Koala
{

    PlanarSSSP::PlanarSSSP(NetworKit::Graph &graph, NetworKit::node source,
                           NetworKit::node target)
        : graph(graph), source(source), target(target) {}

    NetworKit::count PlanarSSSP::getSSSDistanceToTarget()
    {
        assureFinished();
        return distanceToTarget;
    }

} /* namespace Koala */
