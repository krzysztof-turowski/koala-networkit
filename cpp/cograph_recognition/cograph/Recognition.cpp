#include <recognition/CographRecognition.hpp>
#include <list>
#include <set>

#include <graph/GraphTools.hpp>
namespace Koala {
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
    vector<CoNode>outEdges;//neighbours of cur vertex in G
    public:
    CoNode(Type type, int number): type(type), number(number),marked(Marked::UNMARKED),md(0),d(0),head_of_list_of_children(nullptr),
    next(nullptr),prev(nullptr), parent(nullptr),outEdges(vector<CoNode>(0)){

    }
    void addchild(CoNode *x){
        if(head_of_list_of_children == nullptr){
            head_of_list_of_children = x;
            x -> prev = x -> next = nullptr;
        } else{
            head_of_list_of_children -> prev = x;
            x -> next = head_of_list_of_children;
            x -> prev = nullptr;
            head_of_list_of_children = x;
        }
        x -> parent = this;
        d++;
    }
    void setoutEdges(vector<CoNode>outEdges){
        this.outEdges = outEdges;
    }
    vector<CoNode> getoutEdges(){
        return outEdges;
    }
    CoNode* getParent(){
        return parent;
    }
    CoNode* setParent(CoNode *parent){
    this.parent = parent;
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
    md = 0;
    }
    void mark(){
    marked = Marked::MARKED;}
    void unmark(){
        marked = Marked::MARKED_AND_UNMARKED;
    }
}
class CoTree{
    private:
    CoNode* root;
    public:
    CoTree(CoNode* root):root(root){
    }
    CoNode* getRoot(){
        return root;
    }
    void setRoot(CoNode* R){
        root = R;
    }

};
CoTree T;
NetworKit::Graph G;
int mark_count = 0;
int mark_and_unmarked_count = 0;
int mark_ever_count = 0;
void unmark(queue<CoNode> &marked_with_d_equal_to_md){
    CoNode u = marked_with_d_equal_to_md.front();
    marked_with_d_equal_to_md.pop();
    u.unmark();
    mark_count--;
    mark_and_unmarked_count++;
    u.set_md(0);
    if(u != T.getRoot()){
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
void Mark( CoNode &x){
    queue<CoNode> marked_with_d_equal_to_md;
    mark_count = 0;
    mark_and_unmarked_count = 0;
    mark_ever_count = 0;
    for(auto u : T.outEdges(x)){
        u.mark();
        mark_ever_count++;
        mark_count++;
        marked_with_d_equal_to_md.push(u);
    }
    while(!marked_with_d_equal_to_md.empty()){
        unmark( marked_with_d_equal_to_md);
    }
        if(mark_count){
            if(T.getRoot().d == 1){
                T.getRoot().mark();
            }
        }

}

void Reset_All_CoNodes(CoNode *x){
    x -> unmark_for_new_iteration();
    auto y = x -> head_of_list_of_children;
    while(y != nullptr){
        Reset_All_CoNodes(y);
        y = y -> next;
    }
}
CoNode Find_Lowest( ){

}
void Insert_x_to_CoTree(CoNode *u, CoNode *x){
    if(u -> )
}
void Cograph_Recognition(NetworKit::Graph &graph){
    G = graph;
    vector<node>vertex;
    vector<CoNode>covertex;
    vector<int>pos(G.z + 1);
    for(int i = 0; i <= G.z; i++){
        if(G.exists(i)){
            pos[i] = vertex.size();
            vertex.push_back(i);
            CoNode C = CoNode(Type::VERTEX,i);
            covertex.push_back(C);
        }
    }

    for(int i = 0; i <= G.z; i++){
            if(G.exists(i)){
                vector<CoNode> vec;
                for(auto u : G.outEdges(i)){
                    vec.push_back(covertex[pos[u]]);
                }
                covertex[pos[i]].setoutEdges(vec);
            }
        }
    CoNode R(Type::ZEROONE, 1);
    T = CoTree(&R);
    if(vertex.size() == 0){
        return;
    }
    if(vertex.size() == 1){
        R.addchild(&covertex[0]);
        return;
    }
    if(G.hasEdge(vertex[0], vertex[1])){
        R.addchild(&covertex[0]);
        R.addchild(&covertex[1]);
    } else{
        CoNode N(Type::ZEROONE, 0);
        R.addchild(&N);
        N.addchild(&covertex[0]);
        N.addchild(&covertex[1]);
    }
    for(int i = 2; i < vertex.size()){
        Reset_All_CoNodes(&R);
        Mark(covertex[i]);
        if(R.marked == Marked::MARKED_AND_UNMARKED){//all nodes of T were marked and unmarked <=> R is marked and unmarked
            R.addchild(&covertex[i]);
        } else if(mark_ever_count == 0){
            if(R.d == 1){
                (R.head_of_list_of_children) -> addchild(&covertex[i]);
            } else{
                CoNode R1(Type::ZEROONE, 1);
                CoNode R2(Type::ZEROONE, 0);
                R1.addchild(&R2);
                R2.addchild(&R);
                R2.addchild(&x);
            }
        } else{
            CoNode u = Find_Lowest();
            Insert_x_to_CoTree(&u, &covertex[i]);
        }
    }
}


}