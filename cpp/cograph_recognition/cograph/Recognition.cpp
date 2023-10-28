#include <recognition/CographRecognition.hpp>
#include <list>
#include <set>

#include <graph/GraphTools.hpp>
namespace Koala {
/*class Parameters_of_one_iteration_of_algorithm{
    private:
    vector<bool>Marked;
    vector<int> md, d;
    vector<vector<NetworKit::node>> marked_and_unmarked;
    int size;
    queue<int> marked_with_d_equal_to_md;
    public:
    Parameters_of_one_iteration_of_algorithm(int n): size(n){
        Marked.resize(n);
        md.resize(n);
        d.resize(n);
        marked_and_unmarked.resize(n);
        marked_with_d_equal_to_md.clean();
    }

};*/
enum class Type{
ZEROONE,
VERTEX
};
enum class Marked{
    UNMARKED,
    MARKED,
    MARKED_AND_UNMARKED
}
class CoNode{
    private:
    Type type;
    int number;
    Marked marked;
    int md, d;
    //md is the current number of children, which have been both "marked" and "unmarked"
    //d is the current number of children

    CoNode* head_of_list_of_children;
    CoNode *next, *prev;//in list of children of its parent
    CoNode *parent;
    public:
    CoNode(Type type, int number, CoNode *parent): type(type), number(number),marked(Marked::UNMARKED),md(0),d(0),head_of_list_of_children(nullptr),
    next(nullptr),prev(nullptr), parent(parent){

    }
    CoNode* getParent(){
        return parent;
    }
    void set_md(int md){
        this.md = md;
    }
    void inc_md(){
        this.md++;
    }
    int get_md(){
        return md;
    }
    void set_d(int d){
            this.d = d;
        }
        void inc_d(){
            this.d++;
        }
        int get_d(){
            return d;
        }
    Marked Marked_or_not(){
        return marked;
    }
    void unmark_for_new_iteration(){
    marked = Marked::UNMARKED;
    }
    void mark(){
    marked = Marked::MARKED;}
    void unmark(){
        marked = Marked::MARKED_AND_UNMARKED;
    }
}
class CoTree{
    private:
    CoNode root;
    public:
    CoTree(CoNode root):root(root){
    }
    CoNode getRoot(){
        return root;
    }

};
//Parameters_of_one_iteration_of_algorithm param;
CoTree T;
NetworKit::Graph G;
void unmark(queue<CoNode> &marked_with_d_equal_to_md){
    CoNode R = T.getRoot();
    CoNode u = marked_with_d_equal_to_md.front();
    marked_with_d_equal_to_md.pop();
    u.set_md(0);
    if(u != R){
        w = u.getParent();
        w.inc_md();
        if(w.get_md() == w.get_d()){
            marked_with_d_equal_to_md.push(w);
            auto l = u -> next;
            if(u.prev != nullptr){
                (u.prev) -> next = l;
                u.prev = nullptr;
                u.next = w.head_of_list_of_children;
                w.head_of_list_of_children.prev = u;
                w.head_of_list_of_children = u;
            }
            //else u is head
         }
    }
}
void Mark( NetworKit::node &x){
    CoNode R = T.getRoot();
    queue<CoNode> marked_with_d_equal_to_md;
    for(auto u : G.outEdges){//
        u.mark();
        marked_with_d_equal_to_md.push(u);
    }
    while(!marked_with_d_equal_to_md.empty()){
        unmark( marked_with_d_equal_to_md);
    }
    auto nodes = NodeRange(G);
    for(auto u : nodes){
        if(Marked[u]){
            if(d[R] == 1){
                Marked[R] = 1;
            }
            break;
        }
    }
}
void Find_Lowest( NetworKit::node &R){

}
void Cograph_Recognition(NetworKit::Graph &graph){
    G = graph;
    /*for(auto u : NodeRange(graph)){
        Mark( u,R);
    }*/
}


}