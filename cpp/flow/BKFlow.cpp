#include <climits>
#include <algorithm>
#include <queue>
#include <flow/MaximumFlow.hpp>

using edge = std::pair<NetworKit::node, NetworKit::node>;
#define SOURCE 0
#define TARGET 1
#define FREE 2
namespace Koala {

const NetworKit::node NO_PARENT = std::numeric_limits<NetworKit::node>::max();

edge BKFlow::rev(const edge &p) {
    return std::make_pair(p.second, p.first);
}

void BKFlow::initialize(){
    V = graph->numberOfNodes();
    graph->forNodes([&](NetworKit::node v) {
        tree[v]=FREE;
        parent[v]=NO_PARENT;
    });
    tree[source]=SOURCE;
    tree[target]=TARGET;
    active.push(source);
    active.push(target);
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

int  BKFlow::tree_capacity(NetworKit::node p, NetworKit::node q){
    edge forward = std::make_pair(p, q);
    edge backward = std::make_pair(q, p);
    if(tree[p] == SOURCE){
        return capacity[forward] - flow[forward];
    }
    if(tree[p] == TARGET){
        return capacity[backward] - flow[backward];
    }
    return 0;
}

bool  BKFlow::grow(){
    
    while(!active.empty()){
        NetworKit::node v = active.front();
        active.pop();
        bool foundPath = false;
        
        if (tree[v] == FREE) continue; // in case invalid node stayed on the queue
        graph->forNeighborsOf(v, [&](NetworKit::node w) {
            if (foundPath) return;
            edge e = std::make_pair(v,w);
            if(tree_capacity(v,w) > 0){
                if(tree[w] == FREE){
                    tree[w] = tree[v];
                    parent[w] = v;
                    active.push(w);
                }
                else if(tree[w] != tree[v]){
                    spath = tree[v]==SOURCE ? v : w;
                    tpath = tree[w]==TARGET ? w : v;
                    foundPath=true;
                    active.push(v);
                    return;
                }
            }
        });
        if (foundPath) return true;
    }
    return false;
}

int BKFlow::augment() {
    std::vector<edge> path;
    int bottleneck = INT_MAX;

    NetworKit::node u = spath;
    while (u != source) {
        NetworKit::node p = parent[u];
        edge e = std::make_pair(p, u);
        bottleneck = std::min(bottleneck, capacity[e] - flow[e]);
        path.push_back(e);
        u = p;
    }

    edge middle = std::make_pair(spath, tpath);
    bottleneck = std::min(bottleneck, capacity[middle] - flow[middle]);
    path.push_back(middle);

    u = tpath;
    while (u != target) {
        NetworKit::node p = parent[u];
        edge e = std::make_pair(u, p);  // reverse direction
        bottleneck = std::min(bottleneck, capacity[e] - flow[e]);
        path.push_back(e);
        u = p;
    }

    for (edge e : path) {
        flow[e] += bottleneck;
        flow[rev(e)] -= bottleneck;
    }

    u = spath;
    while (u != source) {
        NetworKit::node p = parent[u];
        edge e = std::make_pair(p, u);
        if (capacity[e] - flow[e] == 0) {
            if (tree[p] == SOURCE && tree[u] == SOURCE) {
                parent[u] = NO_PARENT;
                orphan.push(u);
            }
        }
        u = p;
    }

    u = tpath;
    while (u != target) {
        NetworKit::node p = parent[u];
        edge e = std::make_pair(u, p); // reverse direction
        if (capacity[e] - flow[e] == 0) {
            if (tree[p] == TARGET && tree[u] == TARGET) {
                parent[u] = NO_PARENT;
                orphan.push(u);
            }
        }
        u = p;
    }
    return bottleneck;
}

bool BKFlow::origin(NetworKit::node v){
    NetworKit::node u = v;
    while(true){
        if(parent[u]==source || parent[u] == target)return true;
        if(parent[u]==NO_PARENT)return false;
        u = parent[u];
    }
}

void BKFlow::adopt() {
    while (!orphan.empty()) {
        
        NetworKit::node p = orphan.front();
        orphan.pop();
        bool found_new_parent = false;

        graph->forNeighborsOf(p, [&](NetworKit::node q) {
            if(found_new_parent)return;

            if (tree[q] != tree[p]) return;
            if (tree_capacity(q, p) <= 0) return;
            if (!origin(q)) return;

            parent[p] = q;
            found_new_parent = true;
        });

        if (!found_new_parent) {

            graph->forNeighborsOf(p, [&](NetworKit::node q) {
                if (tree[q] != tree[p]) return;

                if (tree_capacity(q, p) > 0) {
                    active.push(q); 
                }

                if (parent[q] == p) {
                    parent[q] = NO_PARENT;
                    orphan.push(q);
                }
            });
            tree[p] = FREE;
        }
        
    }
}

void BKFlow::run() {

    initialize();
    int totalflow = 0;
    int iteration = 0;

    while (true) {

        if (!grow()) {
            break;
        }

        int f = augment();
        totalflow += f;
        adopt();
  
    }

    flow_size = totalflow;
    hasRun = true;

}

}