#include <climits>
#include <algorithm>
#include <queue>
#include <flow/MaximumFlow.hpp>

using edge = std::pair<NetworKit::node, NetworKit::node>;

namespace Koala {
    

const int UNREACHABLE = -2;
    
edge MKMFlow::rev(const edge &p) {
    return std::make_pair(p.second, p.first);
}

void MKMFlow::initialize(){

    V = graph->numberOfNodes();
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        
        auto p = std::make_pair(u, v);
        if(graph->addEdge(v,u,0,true)){
            capacity[rev(p)]=0;
        }
        flow[p] = 0;
        flow[rev(p)]=0;
        capacity[p]=w;
    });
}

bool MKMFlow::buildLevelGraph(){

    graph->forNodes([&](NetworKit::node v) {
        level[v]=UNREACHABLE;
    });

    std::queue<NetworKit::node> q;
    q.push(source);
    level[source] = 0;

    while (!q.empty()) {

        NetworKit::node u = q.front();
        q.pop();

        graph->forNeighborsOf(u, [&](NetworKit::node w) {

            auto e = std::make_pair(u, w);
            if (level[w] == UNREACHABLE && flow[e] < capacity[e]) {
                level[w] = level[u] + 1;
                q.push(w);
            }
        });

    }
    graph_stage = NetworKit::Graph(*graph);
    return level[target] != UNREACHABLE;

}

void MKMFlow::computePotential(){

    graph->forNodes([&](NetworKit::node v) {

        inPotential[v]=0;
        outPotential[v]=0;

    });

    graph->forEdges([&](NetworKit::node x, NetworKit::node y, NetworKit::edgeweight weight){
        edge e = std::make_pair(x,y);
        if(level[e.first] + 1 == level[e.second]){

            if (capacity[e] > flow[e]) {
                outPotential[e.first] += capacity[e] - flow[e];
                inPotential[e.second] += capacity[e] - flow[e];
            }
        }

    });
           
    inPotential[source]=INT_MAX;
    outPotential[target]=INT_MAX;

}

void MKMFlow::pushForward(NetworKit::node u, int f){
    if (u == target){ 
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node,int>to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v]=0;
    });

    to_push[u] += f;
    q.push(u);

    while(!q.empty()){

        NetworKit::node v = q.front();
        q.pop();
        if(to_push[v]==0){
            continue;
        }
        std::vector<edge> edgesToRemove;
        graph_stage.forEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight){ 
            edge e = std::make_pair(v,w);
            if(level[v] + 1 != level[e.second]){
                return;
            }
            long long can_be_pushed = std::min(capacity[e]-flow[e],to_push[v]);
            if(can_be_pushed==0){
                return;
            }
            if(e.second != target){
                q.push(e.second);
            }
            flow[e]+=can_be_pushed;
            flow[rev(e)]-=can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edgesToRemove.push_back(e);
            }
            inPotential[e.second]-=can_be_pushed;
            outPotential[v]-=can_be_pushed;
            to_push[v]-=can_be_pushed;
            to_push[e.second]+=can_be_pushed;
        });
        for (const edge& e : edgesToRemove) {
            graph_stage.removeEdge(e.first, e.second);
            graph_stage.removeEdge(e.second, e.first); 
        }
    }
}

void MKMFlow::pushBackward(NetworKit::node u, int f){
    if (u == source){ 
        return;
    }
    std::queue<NetworKit::node> q;
    std::map<NetworKit::node,int>to_push;

    graph_stage.forNodes([&](NetworKit::node v) {
        to_push[v]=0;
    });

    to_push[u] += f;
    q.push(u);

    while(!q.empty()){

        NetworKit::node v = q.front();
        q.pop();
        if(to_push[v]==0){
            continue;
        }
        std::vector<edge> edgesToRemove;
        graph_stage.forInEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight){ 
            edge e = std::make_pair(w,v);
            if(level[v] - 1 != level[e.first]){
                return;
            }
            long long can_be_pushed = std::min(capacity[e]-flow[e],to_push[v]);
            if(can_be_pushed==0){
                return;
            }
            if(e.first != source){
                q.push(e.first);
            }
            flow[e]+=can_be_pushed;
            flow[rev(e)]-=can_be_pushed;
            if (capacity[e] - flow[e] == 0) {
                edgesToRemove.push_back(e);
            }
            outPotential[e.first]-=can_be_pushed;
            inPotential[v]-=can_be_pushed;
            to_push[v]-=can_be_pushed;
            to_push[e.first]+=can_be_pushed;
        });
        for (const edge& e : edgesToRemove) {
            graph_stage.removeEdge(e.first, e.second);
            graph_stage.removeEdge(e.second, e.first); 
        }
    }
}

void MKMFlow::deleteNode(NetworKit::node v){

     graph_stage.forInEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight){
        edge e = std::make_pair(w,v);
        if(level[v] - 1 == level[e.first]){
            if (capacity[e] > flow[e]) {

                outPotential[e.first] -= capacity[e] - flow[e];
            }
        }
    });
    graph_stage.forEdgesOf(v, [&](NetworKit::node w, NetworKit::edgeweight weight){
        edge e = std::make_pair(v,w);
        if(level[v] + 1 == level[e.second]){
            if (capacity[e] > flow[e]) {

                inPotential[e.second] -= capacity[e] -flow[e];
            }
        }
    });
    level[v]=UNREACHABLE;
    graph_stage.removeNode(v);
}

void MKMFlow::run(){
    initialize();
    int totalFlow = 0;
    while (buildLevelGraph()) {

        computePotential();
        while(true){

            NetworKit::node u;
            int minimum=INT_MAX;

            graph_stage.forNodes([&](NetworKit::node v) {

                if(level[v]==UNREACHABLE){
                    return;
                }
                int pot = std::min(outPotential[v],inPotential[v]);
                if(pot < minimum){
                    u = v;
                    minimum = pot;
                }
            });
            if(minimum==0){
                deleteNode(u);
                continue;
            }
            if(minimum==INT_MAX){
                break;
            }

            pushForward(u,minimum);
            pushBackward(u,minimum);
            totalFlow+=minimum;
            deleteNode(u);
        }

    }

    flow_size = totalFlow;
    hasRun = true;
}

}