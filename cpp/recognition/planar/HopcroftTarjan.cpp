/*
 * HopcroftTarjan.cpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <list>
#include <set>
#include <iostream>
#include <recognition/planar/HopcroftTarjan.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace Koala {

    HopcroftTarjan::HopcroftTarjan(NetworKit::Graph &graph, bool embedding): IPlanarGraphRecognition(graph, embedding){ }

    void HopcroftTarjan::dfs(NetworKit::Graph &g,
                            NetworKit::node v, 
                            std::vector<bool> & reached
                            )
    {
        reached[v] = true;
        for(auto u: g.neighborRange(v)){
            if(!reached[u]){
                dfs(g, u, reached);
            }
        }
    }

    void HopcroftTarjan::dfs_in_make_biconnected_graph(NetworKit::Graph &g, 
                                                    NetworKit::node v, 
                                                    int & dfs_count, 
                                                    std::vector<bool> & reached,
                                                    std::vector<NetworKit::node> & dfsnum,
                                                    std::vector<NetworKit::node> &lowpt, 
                                                    std::vector<NetworKit::node> &parent)
    {
        NetworKit::node w = INT_MAX;
        dfsnum[v] = dfs_count++;
        lowpt[v] = dfsnum[v];
        reached[v] = true;
        std::vector<int> neighbors;
        for(auto u: g.neighborRange(v)) {
            neighbors.push_back(u);
        }
        for(auto u: neighbors){
            if(!reached[u]){
                if(w == INT_MAX){
                    w = u;
                }
                parent[u] = v;
                dfs_in_make_biconnected_graph(g, u, dfs_count, reached, dfsnum, lowpt, parent);
                if(lowpt[u] == dfsnum[v]){
                    if(u == w && parent[v] != INT_MAX) {
                        g.addEdge(u, parent[v]);
                        g.addEdge(parent[v], u);
                    }
                    if(u != w){
                        g.addEdge(u, w);
                        g.addEdge(w, u);
                    }
                }
                lowpt[v] = std::min(lowpt[v], lowpt[u]);
            } else {
                lowpt[v] = std::min(lowpt[v], dfsnum[u]);
            }
        }
    }

    void HopcroftTarjan::make_biconnected_graph(NetworKit::Graph &g) {
        // make it connected
        NetworKit::node u = 0;
        std::vector<bool> reached(N);
        for(auto v: g.nodeRange()){
            if(!reached[v]) {
                dfs(g, v, reached);
                if(u != v) {
                   g.addEdge(u, v);
                   g.addEdge(v, u);
                }
            }
        }
        // make it biconnected
        fill(reached.begin(), reached.end(), false);
        std::vector<NetworKit::node> lowpt(N);
        std::vector<NetworKit::node> parent(N, INT_MAX);
        std::vector<NetworKit::node> dfsnum(N);
        int dfs_count = 0;

        dfs_in_make_biconnected_graph(g, u, dfs_count, reached, dfsnum, lowpt, parent);
    }

    void HopcroftTarjan::dfs_in_reorder(NetworKit::Graph &g, 
                                        NetworKit::node v, 
                                        std::vector<NetworKit::Edge> & Add,
                                        int & dfs_count, 
                                        std::vector<bool> & reached,
                                        std::vector<NetworKit::node> & dfsnum,
                                        std::vector<NetworKit::node> &lowpt, 
                                        std::vector<NetworKit::node> &lowpt2, 
                                        std::vector<NetworKit::node> &parent)
    {
        dfsnum[v] = dfs_count++;
        lowpt2[v] = lowpt[v] = dfsnum[v];
        reached[v] = true;
        for(auto w: g.neighborRange(v)){
            if(!reached[w]){
                Add.push_back(NetworKit::Edge(v, w));
                parent[w] = v;
                dfs_in_reorder(g, w, Add, dfs_count, reached, dfsnum, lowpt, lowpt2, parent);
                lowpt[v] = std::min(lowpt[v], lowpt[w]);
            } else {
                lowpt[v] = std::min(lowpt[v], dfsnum[w]);
                if(dfsnum[w] < dfsnum[v] && parent[v] != w){
                    Add.push_back({v, w});
                }
            }
        }
        for(auto w: g.neighborRange(v)){
            if(parent[w] == v){
                if(lowpt[w] != lowpt[v]) {
                    lowpt2[v] = std::min(lowpt2[v], lowpt[w]);
                }
                lowpt2[v] = std::min(lowpt2[v], lowpt2[w]);
            } else {
                if(lowpt[v] != dfsnum[w]){
                     lowpt2[v] = std::min(lowpt2[v], dfsnum[w]);
                }
            }
        }
    }


    void HopcroftTarjan::reorder(NetworKit::Graph &g, 
                std::vector<NetworKit::node> & parent, 
                std::vector<NetworKit::node> & dfsnum,
                std::vector<NetworKit::Edge> & order
                )
    {
        NetworKit::node v = 0;
        std::vector<bool> reached(N);
        int dfs_count = 0;
        std::vector<NetworKit::Edge> Add;
        std::vector<NetworKit::node> lowpt(N), lowpt2(N);
        dfs_in_reorder(g, v, Add, dfs_count, reached, dfsnum, lowpt, lowpt2, parent);
        std::vector<NetworKit::Edge> buckets[N*2+5];
        for(auto u: Add){
            int cost = ((dfsnum[u.v] < dfsnum[u.u]) ? 2 * dfsnum[u.v] :
            ((lowpt2[u.v] >= dfsnum[u.u]) ? 2 * lowpt[u.v] : 2 * lowpt[u.v] + 1));
            buckets[cost].push_back(u);
        }
        g = NetworKit::Graph(N, true, true);
        order = std::vector<NetworKit::Edge>(Add.size());
        int id = 0;
        for(int i = 0; i <= 2*N; i++){
            for(auto u: buckets[i]){
                g.addEdge(u.u, u.v, id);
                order[id] = NetworKit::Edge(u.u, u.v);
                id++;
            }
        }
    }

    bool HopcroftTarjan::strong_planar(NetworKit::Graph &g, 
                       NetworKit::Edge e, 
                       std::list<int> & Att,
                       std::vector<int> & alpha,
                       std::vector<NetworKit::node> & dfsnum,
                       std::vector<NetworKit::node> & parent
                    )
    {
        NetworKit::node  y = e.v, t = *g.neighborRange(y).begin(), x = e.u, wk = y;
        while(dfsnum[t] > dfsnum[wk]){
            wk = t;
            t = *g.neighborRange(wk).begin();
        }
        NetworKit::node w0 = t;
        NetworKit::node w = wk;
        std::stack<Block*> S;
        while(w != x) {
            int count = 0;
            for(auto u : g.weightNeighborRange(w)) {
                count++;
                if(count != 1) {
                    std::list<int> A;
                    if(dfsnum[w] < dfsnum[u.first]) {
                        if(!strong_planar(g, NetworKit::Edge(w, u.first), A, alpha, dfsnum, parent)){
                            while(!S.empty()){
                                delete S.top();
                                S.pop();
                            } 
                            return false;
                        }
                    } else {
                        A.push_back(dfsnum[u.first]);
                    }

                    Block * b = new Block(NetworKit::WeightedEdge(w, u.first, u.second), A);
                    while(true){
                        if(b->left_interlace(S)){
                            (S.top())->flip();
                        }
                        if(b->left_interlace(S)){
                            delete b;
                            while(!S.empty()){
                                delete S.top();
                                S.pop();
                            }
                            return false;
                        }
                        if(b->right_interlace(S)){
                            b->combine(S.top());
                            S.pop();
                        } else {
                            break;
                        }
                    }
                    S.push(b);
                }
            }
            while(!S.empty() && (S.top())->clean(dfsnum[parent[w]], alpha, dfsnum)) {
                delete S.top();
                S.pop();
            }
            w = parent[w];
        } 
        Att.clear();
        while(!S.empty()){
            Block *B = S.top(); S.pop();
            if(!(B->empty_Latt()) && !(B->empty_Ratt()) && (B->head_of_Latt() > dfsnum[w0]) && (B->head_of_Ratt() > dfsnum[w0])) {
                delete B;
                while(!S.empty()) {
                    delete S.top();  S.pop();
                }
                return false;
            }
            B->add_to_Att(Att, dfsnum[w0], alpha);
            delete B;
        }
        if(w0 != x){
            Att.push_back(dfsnum[w0]);
        }

        return true;
    }

    bool HopcroftTarjan::planar() {
        N = graph.numberOfNodes();
        if(N <= 3) {
            return true;
        }
        if(graph.numberOfEdges() > 6 * N - 12) {
            return false;
        }

        // make G copy of Graph and make it biconnected

        NetworKit::Graph g = NetworKit::Graph(N, false, true);
        for(auto i: graph.edgeRange()){
            g.addEdge(i.u, i.v);
            g.addEdge(i.v, i.u);
        }
        make_biconnected_graph(g);
         if(g.numberOfEdges() > 6 * N - 12) {
            return false;
        }
       
       

        //test planarity
        std::vector<NetworKit::node> parent(N, INT_MAX);
        std::vector<NetworKit::node> dfsnum(N);
        std::vector<NetworKit::Edge> order;
        reorder(g, parent, dfsnum, order);
    
        std::vector<int> alpha(g.numberOfEdges());
        {
            std::list<int> Att;
            NetworKit::Edge e(0, *g.neighborRange(0).begin());
            alpha[0] = left;

            if(!strong_planar(g, e, Att, alpha, dfsnum, parent)){
                return false;
            }
        }
        if(embedding) {
            // TODO embedding 
        }
        return true;
       
    }

    void HopcroftTarjan::run() {
        hasRun = true;
        
        if(this->planar()) {
            this->is_planar = State::PLANAR;
        }
    }

}  // namespace Koala
