#include <recognition/CographRecognition.hpp>
#include <list>
#include <set>

#include <graph/GraphTools.hpp>
namespace Koala {
void unmark(NetworKit::Graph &graph, NetworKit::node &u, vector<bool> &Marked, vector<int> &md, vector<int>&d,NetworKit::node &R,queue<int> &marked_with_d_equal_to_md, vector<vector<NetworKit::node>> &marked_and_unmarked, DirectedTree &T){
    marked_with_d_equal_to_md.pop(u);
    md[u] = 0;
    if(u != R){
        w = T.getParent(u);
        md[w] ++;
        if(md[w] == d[w]){
            marked_with_d_equal_to_md.push(w);
            //TODO: insert u at the head of a linked list of marked and unmarked children of w
         }
    }
}
void Mark(NetworKit::Graph &graph, NetworKit::node &x, vector<bool> &Marked, vector<int> &md, vector<int>&d, NetworKit::node &R,vector<vector<NetworKit::node>> &marked_and_unmarked, DirectedTree &T){
    queue<int> marked_with_d_equal_to_md;
    for(auto u : graph.outEdges){
        Marked[u] = true;
        marked_with_d_equal_to_md.push(u);
    }
    while(!marked_with_d_equal_to_md.empty()){
        unmark(graph, u, Marked, md, d, R, marked_with_d_equal_to_md,marked_and_unmarked, T);
    }
    auto nodes = NodeRange(graph);
    for(auto u : nodes){
        if(Marked[u]){
            if(d[R] == 1){
                Marked[R] = 1;
            }
            break;
        }
    }
}
void Find_Lowest(vector<bool> &Marked,vector<int> &md, vector<int>&d, NetworKit::node &R, vector<vector<NetworKit::node>> &marked_and_unmarked){

}
void Cograph_Recognition(NetworKit::Graph &graph){
    vector<int> md, d;
    vector<vector<NetworKit::node>> marked_and_unmarked;
    vector<bool>Marked;
    NetworKit::Graph G = NetworKit::Graph();
    DirectedTree T = DirectedTree(G);
    ///////////
    for(auto u : NodeRange(graph)){
        Mark(graph, u, Marked, md, d, R, marked_and_unmarked, T);
    }
}


}