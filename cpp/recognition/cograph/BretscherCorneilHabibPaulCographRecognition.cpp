/*
 * BretscherCorneilHabibPaulCographRecognition.cpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
// Copyright 2024 milana
#include <graph/GraphTools.hpp>

#include "recognition/CographRecognition.hpp"
#include "recognition/CoTree.hpp"

namespace Koala {

bool BretscherCorneilHabibPaulCographRecognition::isCograph() const {
    assureFinished();
    return is_cograph;
}

void BretscherCorneilHabibPaulCographRecognition::run() {
    hasRun = true;
    auto a = LexBfs();
    auto b = LexBfsMinus(graph, a);//graph'
    auto c = LexBfsMinus(graph, b);
    is_cograph = NeighbourhoodSubsetProperty(b, c);
}

std::vector<NetworKit::node> BretscherCorneilHabibPaulCographRecognition::LexBfs(){

}

std::vector<NetworKit::node> BretscherCorneilHabibPaulCographRecognition::
LexBfsMinus(NetworKit::Graph G, std::vector<NetworKit::node> a){

}

bool BretscherCorneilHabibPaulCographRecognition::NeighbourhoodSubsetProperty(std::vector<NetworKit::node> a, std::vector<NetworKit::node> b){

}

}  // namespace Koala