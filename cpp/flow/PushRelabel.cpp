#include <climits>
#include <algorithm>
#include <flow/MaximumFlow.hpp>

using edge = std::pair<NetworKit::node, NetworKit::node>;

namespace Koala {

edge PushRelabel::rev(const edge &p) {
    return std::make_pair(p.second, p.first);
}
void PushRelabel::initialize(){
    V = graph->numberOfNodes();
    graph->forNodes([&](NetworKit::node v) {
        height[v]=0;
        nextedge[v]=0;
        excess[v]=0;
    });
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        
        auto p = std::make_pair(u, v);
        if(graph->addEdge(v,u,0,true)){
            capacity[rev(p)]=0;
        }
        flow[p] = 0;
        flow[rev(p)]=0;
        capacity[p]=w;
    });
    height[source] = V;
    excess[source] = INT_MAX;
}

void PushRelabel::push(const NetworKit::node& vertex, const edge& e){

    int fl = std::min(capacity[e] - flow[e], excess[vertex]);
    excess[vertex] -= fl;
    excess[e.second] += fl;
    flow[e]+=fl;
    flow[rev(e)]-=fl;

    if(excess[e.second] > 0 && e.second != source && e.second != target){
        q.push(e.second);
    }

}

void PushRelabel::relabel(const NetworKit::node& vertex){

    int minimum = INT_MAX;

    graph->forNeighborsOf(vertex, [&](NetworKit::node w) {

        auto e = std::make_pair(vertex, w);
        if(capacity[e] - flow[e] > 0) {
            minimum = std::min(minimum, height[e.second]);
        }
    });

    if(minimum != INT_MAX) {
        height[vertex] = minimum + 1;
    }

}

void PushRelabel::discharge(const NetworKit::node& vertex){

    while(excess[vertex] > 0){

        while(nextedge[vertex] < graph->degreeOut(vertex)){

            auto u = graph->getIthNeighbor(vertex,nextedge[vertex]);
            auto e = std::make_pair(vertex,u);

            if(capacity[e] - flow[e] > 0 && height[vertex] == height[e.second] + 1){
                push(vertex,e);
                if(excess[vertex] == 0){
                    return;
                } 
            }
            nextedge[vertex]++;
        }

        relabel(vertex);
        nextedge[vertex] = 0;

        if(height[vertex] >= V) break;
    }
}

void PushRelabel::run(){
    //printf("Push Relabel Running\n");
    initialize();

    graph->forNeighborsOf(source, [&](NetworKit::node w) {
        push(source,std::make_pair(source,w));
    });

    while(!q.empty()){
        auto vertex = q.front();
        q.pop();
        discharge(vertex);
    }

    flow_size = excess[target];
    hasRun = true;
}


}