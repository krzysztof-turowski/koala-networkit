/*
 * PlanarGraphRecognition.cpp
 *
 *  Created on: 06.04.2023
 */

#include <recognition/PlanarGraphRecognition.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

namespace Koala {

PlanarGraphRecognition::PlanarGraphRecognition(NetworKit::Graph &graph)
    : graph(std::make_optional(graph)), is_planar(false) { }

bool PlanarGraphRecognition::isPlanar() const {
    assureFinished();
    return is_planar;
}

void PlanarGraphRecognition::run() {
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G(graph->numberOfNodes());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        boost::add_edge(u, v, G);
    });
    is_planar = boost::boyer_myrvold_planarity_test(G);
    hasRun = true;
}

}  // namespace Koala
