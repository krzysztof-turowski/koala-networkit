/*
 * IPlanarGraphRecognition.cpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <recognition/planar/IPlanarGraphRecognition.hpp>

namespace Koala {
    IPlanarGraphRecognition::IPlanarGraphRecognition(NetworKit::Graph &graph, bool embedding = false): graph(graph), embedding(embedding), is_planar(State::NOT_PLANAR) {}

    bool IPlanarGraphRecognition::isPlanar() const {
        assureFinished();
        return is_planar == State::PLANAR;
    }

    void IPlanarGraphRecognition::getEmbedding() const {
        assureFinished();
        return;
    }

}  // namespace Koala
