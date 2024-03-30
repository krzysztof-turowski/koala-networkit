#include <cograph_recognition/CographRecognition.hpp>
#include <list>
#include <graph/GraphTools.hpp>
namespace Koala {
    enum class Type{
        ZERO_ONE,
        VERTEX
    };
    enum class Marked{
        UNMARKED,
        MARKED,
        MARKED_AND_UNMARKED
    };

    class CoNode{
    public:
        Type type;
        int number;
        Marked marked;
        int md, d;
        bool in_graph;

        //d is the current number of children
        //md is the current number of children, which have been both "marked" and "unmarked"

        CoNode* head_of_list_of_children;
        CoNode *next, *prev;//in list of children of its parent
        CoNode *parent;
        std::vector<CoNode*>out_edges;//neighbours of cur vertex in G

        CoNode(Type type, int number): type(type), number(number),marked(Marked::UNMARKED),md(0),d(0),in_graph(false),head_of_list_of_children(nullptr),
        next(nullptr), prev(nullptr),parent(nullptr){

        }

        void AddChild(CoNode *x){
            if(head_of_list_of_children == nullptr){
                head_of_list_of_children = x;
                x -> prev = nullptr;
                x -> next = nullptr;
            } else{
                head_of_list_of_children -> prev = x;
                x -> next = head_of_list_of_children;
                x -> prev = nullptr;
                head_of_list_of_children = x;
            }
            x -> parent = this;
            d++;
        }

        void UnmarkForNewIteration(){
            marked = Marked::UNMARKED;
            md = 0;
        }
        void mark(){
            marked = Marked::MARKED;
        }
        void unmark(){
            marked = Marked::MARKED_AND_UNMARKED;
        }
        std::vector<CoNode*> RemoveWereMarked(){
            auto u = head_of_list_of_children;
            std::vector<CoNode*>vec;
            while(u != nullptr){
                vec.push_back(u);
                d--;
                head_of_list_of_children = u -> next;
                if(head_of_list_of_children != nullptr){
                    head_of_list_of_children -> prev = nullptr;
                }
                u -> prev = nullptr;
                u -> next = nullptr;
                u = head_of_list_of_children;
                if(u == nullptr || u -> marked != Marked::MARKED_AND_UNMARKED){
                    break;
                }
            }
            return vec;
        }
        void RemoveWereNotMarked(){
            auto u = head_of_list_of_children;
            while(u != nullptr && u -> marked == Marked::MARKED_AND_UNMARKED){
                u = u -> next;
            }
            while(u != nullptr){
                d--;
                auto nxt = u -> next;
                auto prv = u -> prev;
                if(nxt != nullptr){
                    nxt -> prev = prv;
                }
                if(prv != nullptr){
                    prv -> next = nxt;
                }
                u -> prev = nullptr;
                u -> next = nullptr;
                u = nxt;
            }
        }
    };
    class CoTree{
    private:
        std::vector<CoNode*>save;
    public:
        CoNode* root;
        explicit CoTree(CoNode* root):root(root){
            save.push_back(root);
        }
        void Add(CoNode* x){
            save.push_back(x);
        }
        void Clear(){
            for(auto u : save){
                delete u;
            }
        }
    };
    CoTree* T;
    NetworKit::Graph G;
    int mark_count = 0;
    int mark_and_unmarked_count = 0;
    int mark_ever_count = 0;
    void Unmark(std::queue<CoNode*> &marked_with_d_equal_to_md){
        CoNode* u = marked_with_d_equal_to_md.front();
        marked_with_d_equal_to_md.pop();
        u -> unmark();
        mark_count--;
        mark_and_unmarked_count++;
        u -> md = 0;
        if(u != T -> root){
            auto w = u -> parent;
            w -> md++;
            if(w -> marked == Marked::UNMARKED){
                mark_count++;
            }
            w -> mark();
            if(w -> md == w -> d){
                marked_with_d_equal_to_md.push(w);
            }
            auto nxt = u -> next;
            auto prv = u -> prev;
            auto head = w -> head_of_list_of_children;
            if(prv != nullptr){
                prv -> next = nxt;
                if(nxt != nullptr){
                    nxt -> prev = prv;
                }
                u -> prev = nullptr;
                u -> next = head;
                head -> prev = u;
                w -> head_of_list_of_children = u;
            }//else u is head
        }
    }
    void Mark( CoNode *x){
        std::queue<CoNode*> marked_with_d_equal_to_md;
        mark_count = 0;
        mark_and_unmarked_count = 0;
        mark_ever_count = 0;
        for(auto u : x -> out_edges){//!!only neighbours which are already in graph
            if(!(u -> in_graph)){
                continue;
            }
            u -> mark();
            mark_ever_count++;
            mark_count++;
            marked_with_d_equal_to_md.push(u);
        }
        while(!marked_with_d_equal_to_md.empty()){
            Unmark( marked_with_d_equal_to_md);
        }
        if(mark_count){
            if(T -> root -> d == 1){
                T -> root -> mark();
            }
        }

    }

    void ResetAllCoNodes(CoNode *x,int level=0){
        x -> UnmarkForNewIteration();
        CoNode* y = x -> head_of_list_of_children;
        while(y != nullptr){
            ResetAllCoNodes(y,level+1);
            y = y -> next;
        }
    }
    void MadeQueueOfMarked(CoNode *x, std::queue<CoNode*>&q){
        if(x -> marked == Marked::MARKED){
            q.push(x);
        }
        auto y = x -> head_of_list_of_children;
        while(y != nullptr){
            MadeQueueOfMarked(y, q);
            y = y -> next;
        }
    }
    void Rec(CoNode* x, int level=0){
          std::cout << x << "   " << x -> number << " " << level << " ";
         if(x -> type == Type::ZERO_ONE){
             std::cout<<"zeroone"<<std::endl;
         }
         else {
             std::cout<<"normal"<<std::endl;
         }
        auto y = x -> head_of_list_of_children;
        while(y != nullptr){
            Rec(y,level+1);
            y = y -> next;
        }
    }

    CoNode *FindLowest( CographRecognition::State & error){
        auto* y = new CoNode(Type::ZERO_ONE, 2);
        T -> Add(y);
        CoNode *u, *w, *t;
        if(T -> root -> marked == Marked::UNMARKED){
            error = CographRecognition::State::GRANDPARENT_IS_NOT_IN_SET;
            return y;
        }
        if(T -> root -> md != T -> root -> d - 1) {
            y = T -> root;
        }
            T -> root -> unmark();
            T -> root -> md = 0;
            w = T -> root;
            //u = w;
        std::queue<CoNode*>q;
        MadeQueueOfMarked(T -> root, q);
        while(!q.empty()){
            u = q.front();
            q.pop();
            if(u -> marked != Marked::MARKED){
                continue;
            }
            if(y -> number != 2){//1 or 2
                error = CographRecognition::State::CONTAINS_0_NODE;
                return y;
            }
            if(u -> number == 1){
                if(u -> md != u -> d - 1) {
                    y = u;
                }
                if (u -> parent -> marked == Marked::MARKED) {//1 and 6
                    error = CographRecognition::State::CONTAINS_0_NODE;
                    return y;
                } else {
                    t = u -> parent -> parent;
                }
            } else{
                y = u;
                t = u -> parent;
            }
            u -> unmark();
            u -> md = 0;
            while(t != w){
                if(t == T -> root){//4
                    error = CographRecognition::State::NO_ONE_PATH;
                    return y;
                }
                if(t -> marked != Marked::MARKED){//3 or 5 or 6
                    error = CographRecognition::State::GRANDPARENT_IS_NOT_IN_SET;//!
                    return y;
                }
                if(t -> md != t -> d - 1){//2
                    error = CographRecognition::State::EXISTS_1_NODE_NOT_PROPERLY_MARKED;
                    return y;
                }
                if(t -> parent -> marked == Marked::MARKED){//1
                    error = CographRecognition::State::CONTAINS_0_NODE;
                    return y;
                }
                t -> unmark();
                t -> md = 0;
                t = t -> parent -> parent;
            }
            w = u;
        }
        return w;
    }



    std::vector<CoNode *> GetWereMarked(CoNode *u){
        auto x = u -> head_of_list_of_children;
        std::vector<CoNode *> a;
        while(x != nullptr && x -> marked == Marked::MARKED_AND_UNMARKED){
            a.push_back(x);
            x = x -> next;
        }
        return a;
    }


    CoNode* GetLastFromChildren(CoNode* u){
        auto x = u -> head_of_list_of_children;
        while(x != nullptr && x -> marked == Marked::MARKED_AND_UNMARKED){
            x = x -> next;
        }
        return x;
    }

    void InsertXToCoTree(CoNode *u, CoNode *x){
        std::vector<CoNode*>a;
        int u_number = u -> number;
        a = GetWereMarked(u);
        if((a.size() == 1 && u_number == 0) || (u -> d - a.size() == 1 && u_number == 1)){
            CoNode* w = a[0];
            if(u_number == 1){
                w = GetLastFromChildren(u);
            }
            if(w -> type == Type::VERTEX){
                auto* y = new CoNode(Type::ZERO_ONE, u_number ^ 1);
                T -> Add(y);
                if(u_number == 0){
                    u -> RemoveWereMarked();
                }
                else {
                    u -> RemoveWereNotMarked();
                }
                u -> AddChild(y);
                y -> AddChild(x);
                y -> AddChild(w);
            } else{
                w -> AddChild(x);
            }
        } else{
            auto vec = u -> RemoveWereMarked();
            auto *y = new CoNode(Type::ZERO_ONE, u_number);
            T -> Add(y);
            for(auto v : vec){
                y -> AddChild(v);
            }
            if(u_number == 1){
                auto nxt = u -> next;
                auto prv = u -> prev;
                if(prv != nullptr){
                    prv -> next = y;
                }
                if(nxt != nullptr){
                    nxt -> prev = y;
                }
                y -> prev = prv;
                y -> next = nxt;
                if(prv == nullptr && u -> parent != nullptr){
                    u -> parent -> head_of_list_of_children = y;
                }
                if(u -> parent != nullptr){
                    y -> parent = u -> parent;
                }
                else{
                    T -> root = y;
                }
                auto* z = new CoNode(Type::ZERO_ONE, 0);
                T -> Add(z);
                y -> AddChild(z);
                z -> AddChild(x);
                z -> AddChild(u);
            } else{
                auto *z = new CoNode(Type::ZERO_ONE, 1);
                T -> Add(z);
                u -> AddChild(z);
                z -> AddChild(x);
                z -> AddChild(y);
            }
        }
    }


    CographRecognition::State CographRecognition::Cograph_Recognition(NetworKit::Graph &graph){
        auto *R = new CoNode(Type::ZERO_ONE, 1);
        CoTree Tp(R);
        T = &Tp;
        G = graph;
        std::vector<NetworKit::node>vertex;
        std::vector<CoNode*>covertex;
        std::map<NetworKit::node, int> pos;
        int cnt = 0;
        for(auto i : G.nodeRange()){
            vertex.push_back(i);
            pos[i] = cnt++;
            auto *C = new CoNode(Type::VERTEX,int(i));
            T -> Add(C);
            covertex.push_back(C);
        }
        for(auto i : G.nodeRange()){
            std::vector<CoNode*> vec;
            for(auto u : G.neighborRange(i)){
                vec.push_back(covertex[pos[u]]);
            }
            covertex[pos[i]] -> out_edges = vec;
        }

        if(cnt == 0){
            T -> Clear();
            return State::COGRAPH;
        }
        if(cnt == 1){
            R -> AddChild(covertex[0]);
            T -> Clear();
            return State::COGRAPH;
        }
        if(G.hasEdge(vertex[0], vertex[1])){
            R -> AddChild(covertex[0]);
            R -> AddChild(covertex[1]);
        } else {
            auto *N = new CoNode(Type::ZERO_ONE, 0);
            T -> Add(N);
            R -> AddChild(N);
            N -> AddChild(covertex[0]);
            N -> AddChild(covertex[1]);
        }
        covertex[0] -> in_graph = true;
        covertex[1] -> in_graph = true;
        CoNode* root;

        for(int i = 2; i < cnt; i++){
            root = T -> root;
            ResetAllCoNodes(root);
            Mark(covertex[i]);
            if(root -> marked == Marked::MARKED_AND_UNMARKED){//all nodes of T were marked and unmarked <=> R is marked and unmarked
                root -> AddChild(covertex[i]);
            } else if(mark_ever_count == 0){
                if(root -> d == 1){
                    (root -> head_of_list_of_children) -> AddChild(covertex[i]);
                } else{
                    auto *R1 = new CoNode(Type::ZERO_ONE, 1);
                    auto *R2 = new CoNode(Type::ZERO_ONE, 0);
                    T -> Add(R1);
                    T -> Add(R2);
                    R1 -> AddChild(R2);
                    R2 -> AddChild(root);
                    R2 -> AddChild(covertex[i]);
                    T -> root = R1;
                }
            } else{
                CographRecognition::State error = State::COGRAPH;
                CoNode* u = FindLowest(error);
                if(error != State::COGRAPH){
                    T -> Clear();
                    return error;
                }
                InsertXToCoTree(u, covertex[i]);
            }
            covertex[i] -> in_graph = true;
        }
        T -> Clear();
        return State::COGRAPH;
    }
}
