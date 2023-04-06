/*
 * IndependentSet.cpp
 *
 *  Created on: 06.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <independent_set/IndependentSet.hpp>

namespace Koala {

IndependentSet::IndependentSet(
        NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::set<NetworKit::node>& IndependentSet::getIndependentSet() const {
    assureFinished();
    return out_is;
}

}  /* namespace Koala */
