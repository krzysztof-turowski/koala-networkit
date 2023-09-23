/*
 * PlanarGraphIndependentSet.cpp
 *
 *  Created on: 06.04.2023
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <independent_set/PlanarGraphIndependentSet.hpp>

#include <map>
#include <optional>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>

#include <techniques/baker/bakers_technique.hpp>

#include <techniques/baker/BakerKOuterplanar.hpp>
#include <techniques/baker/Bodlaender.hpp>
#include <techniques/baker/ProblemNode.hpp>

namespace Koala {

BakerPlanarGraphIndependentSet::BakerPlanarGraphIndependentSet(
    NetworKit::Graph &graph, double epsilon)
        : PlanarGraphIndependentSet(graph), epsilon(epsilon) { }

void BakerPlanarGraphIndependentSet::run() {
    auto [G, embedding, outer_face] = get_embedding(*graph);
    int result = baker<IndependentSetNode>(G, *graph, embedding, outer_face);
    for (NetworKit::node u = 0; u < result; u++) {
        independentSet.insert(u);
    }
    hasRun = true;
}

BodlaenderPlanarGraphIndependentSet::BodlaenderPlanarGraphIndependentSet(
    NetworKit::Graph &graph, double epsilon)
        : PlanarGraphIndependentSet(graph), epsilon(epsilon) { }

void BodlaenderPlanarGraphIndependentSet::run() {
    auto [G, embedding, outer_face] = get_embedding(*graph);
    int result = bodlaender<IndependentSetNode>(G, *graph, embedding, outer_face);
    for (NetworKit::node u = 0; u < result; u++) {
        independentSet.insert(u);
    }
    hasRun = true;
}

}  /* namespace Koala */
