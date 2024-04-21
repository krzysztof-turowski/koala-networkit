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

#include <list>

namespace Koala {

bool BretscherCorneilHabibPaulCographRecognition::isCograph() const {
    assureFinished();
    return is_cograph;
}

void BretscherCorneilHabibPaulCographRecognition::run() {
    hasRun = true;
    auto a = LexBfs();
    auto b = LexBfsMinus(true, graph, a);//graph'
    auto c = LexBfsMinus(false, graph, b);
    is_cograph = NeighbourhoodSubsetProperty(b, c);
}

std::vector<NetworKit::node> BretscherCorneilHabibPaulCographRecognition::LexBfs(){
    std::vector<NetworKit::node> a(graph.numberOfNodes());
    std::list<std::list<NetworKit::node>> L;
    std::list<NetworKit::node> first;
    for(auto i : graph.nodeRange()) {
        first.push_back(i);
    }
    L.push_back(first);
    int i = 0;
    while(!L.empty()) {
        auto &l = L.front();
        auto &x = l.front();
        l.pop_front();
        a[x] = i++;
        for(auto it = L.begin(); it != L.end(); it++) {
            auto &j = *it;
            std::list<NetworKit::node> P;
            for(auto u = j.begin(); u != j.end(); u++) {
                if(graph.hasEdge(*u, x)) {
                    P.push_back(*u);
                    j.erase(u);
                }
            }
            if(j.empty()){
                for(auto u : P) {
                    j.push_back(u);
                }
            }
            if(!P.empty()) {
                L.insert(it,P);
                it++;
            }
        }
    }
    return a;
}

std::vector<NetworKit::node> BretscherCorneilHabibPaulCographRecognition::
LexBfsMinus(bool is_complement, NetworKit::Graph G, std::vector<NetworKit::node> a){

}

bool BretscherCorneilHabibPaulCographRecognition::NeighbourhoodSubsetProperty(std::vector<NetworKit::node> a, std::vector<NetworKit::node> b){

}

}  // namespace Koala