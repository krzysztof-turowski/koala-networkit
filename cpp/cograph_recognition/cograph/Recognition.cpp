#include <cograph_recognition/CographRecognition.hpp>
#include <list>
#include <set>
#include <bits/stdc++.h>
#include <graph/GraphTools.hpp>
using namespace std;
namespace Koala {
    enum class Type{
        ZEROONE,
        VERTEX
    };
    enum class Marked{
        UNMARKED,
        MARKED,
        MARKED_AND_UNMARKED
    };

    class CoNode{
    private:
        Type type;
        int number;
        Marked marked;
        int md, d;
        bool in_graph;

        //md is the current number of children, which have been both "marked" and "unmarked"
        //d is the current number of children

        CoNode* head_of_list_of_children;
        CoNode *next, *prev;//in list of children of its parent
        CoNode *parent;
        vector<CoNode*>outEdges;//neighbours of cur vertex in G
    public:
        CoNode(Type type, int number): type(type), number(number),marked(Marked::UNMARKED),md(0),d(0),head_of_list_of_children(nullptr),
                                       next(nullptr),prev(nullptr), parent(nullptr),in_graph(false){

        }
        int getnumber(){
            return number;
        }
        Type gettype(){
            return type;
        }
        void setnumber(int number){
            this -> number = number;
        }
        CoNode *get_head_of_list_of_children(){
            return head_of_list_of_children;
        }
        void set_head_of_list_of_children(CoNode * head_of_list_of_children){
            this -> head_of_list_of_children= head_of_list_of_children;
        }
        CoNode* getnext(){
            return next;
        }
        void setnext(CoNode *next){
            this -> next = next;
        }
        CoNode* getprev(){
            return prev;
        }
        void setprev(CoNode *prev){
            this -> prev = prev;
        }
        void addchild(CoNode *x){
            if(head_of_list_of_children == nullptr){
                head_of_list_of_children = x;
                x ->setprev(nullptr);
                x ->setnext(nullptr);
            } else{
                head_of_list_of_children ->setprev(x);
                x->setnext(head_of_list_of_children);
                x->setprev(nullptr);
                head_of_list_of_children = x;
            }
            x->setParent(this);
            d++;
        }
        void setoutEdges(vector<CoNode*>outEdges){
            this->outEdges = outEdges;
        }
        vector<CoNode*> getoutEdges(){
            return outEdges;
        }
        CoNode* getParent(){
            return parent;
        }
        void setParent(CoNode *parent){
            this->parent = parent;
        }
        void set_md(int md){
            this->md = md;
        }
        void inc_md(){
            this->md = this -> md + 1;
        }
        int get_md(){
            return md;
        }
        void set_d(int d){
            this->d = d;
        }
        void inc_d(){
            this->d++;
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
        vector<CoNode*> remove_were_marked(){
            auto u = head_of_list_of_children;
            vector<CoNode*>vec;
            while(u != nullptr){
                auto save = u;
                vec.push_back(u);
                d--;
                head_of_list_of_children = u -> getnext();
                if(head_of_list_of_children != nullptr)head_of_list_of_children->prev = nullptr;
                u = head_of_list_of_children;
                save ->setprev(nullptr);
                save->setprev(nullptr);
                if(u == nullptr || u -> Marked_or_not() != Marked::MARKED_AND_UNMARKED)break;
            }
            return vec;
        }
        void remove_were_not_marked(){
            auto u = head_of_list_of_children;
            while(u != nullptr && u -> Marked_or_not() == Marked::MARKED_AND_UNMARKED){
                u = u -> getnext();
            }
            while(u != nullptr){
                d--;
                auto save = u;
                auto nxt = u -> getnext();
                auto prv = u -> getprev();
                if(nxt != nullptr)nxt ->setprev(prv);
                if(prv != nullptr)prv ->setnext(nxt);
                u = nxt;
                save ->setprev(nullptr);
                save->setprev(nullptr);
            }
        }
        bool is_in_graph(){
            return in_graph;
        }
        void add_to_graph(){
            in_graph = true;
        }
    };
    class CoTree{
    private:
        CoNode* root;
        vector<CoNode*>save;
    public:
        CoTree(CoNode* root):root(root){
        }
        CoNode* getRoot(){
            return root;
        }
        void add(CoNode* x){
            save.push_back(x);
        }
        void setRoot(CoNode* r){
            root = r;
        }
        void clear(){
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
    void unmark(queue<CoNode*> &marked_with_d_equal_to_md){
        CoNode* u = marked_with_d_equal_to_md.front();
        marked_with_d_equal_to_md.pop();
        u -> unmark();
        mark_count--;
        mark_and_unmarked_count++;
        u -> set_md(0);
        if(u != T -> getRoot()){
            auto w = u -> getParent();
            w -> inc_md();
            if(w -> Marked_or_not() == Marked::UNMARKED)mark_count++;
            w -> mark();
            if(w -> get_md() == w -> get_d()){
                marked_with_d_equal_to_md.push(w);
            }
            auto nxt = u -> getnext();
            auto prv = u -> getprev();
            auto head = w -> get_head_of_list_of_children();
            if(prv != nullptr){
                prv -> setnext(nxt);
                if(nxt != nullptr)nxt ->setprev(prv);
                u -> setprev(nullptr);
                u -> setnext(head);
                head -> setprev (u);
                w -> set_head_of_list_of_children(u);
            }//else u is head
        }
    }
    void Mark( CoNode *x){
        queue<CoNode*> marked_with_d_equal_to_md;
        mark_count = 0;
        mark_and_unmarked_count = 0;
        mark_ever_count = 0;
        for(auto u : x -> getoutEdges()){//!!only neigbours which are already in graph
        if(!(u -> is_in_graph()))continue;
            u -> mark();
            mark_ever_count++;
            mark_count++;
            marked_with_d_equal_to_md.push(u);
        }
        while(!marked_with_d_equal_to_md.empty()){
            unmark( marked_with_d_equal_to_md);
        }
        if(mark_count){
            if(T -> getRoot() -> get_d() == 1){
                T -> getRoot() -> mark();
            }
        }

    }

    void Reset_All_CoNodes(CoNode *x,int level=0){
        x -> unmark_for_new_iteration();
        CoNode* y = x -> get_head_of_list_of_children();
        while(y != nullptr){
            Reset_All_CoNodes(y,level+1);
            y = y -> getnext();
        }
    }
    void Made_Queue_Of_Marked(CoNode *x, queue<CoNode*>&q){
        if(x -> Marked_or_not() == Marked::MARKED){
            q.push(x);
        }
        auto y = x -> get_head_of_list_of_children();
        while(y != nullptr){
            Made_Queue_Of_Marked(y, q);
            y = y -> getnext();
        }
    }
    void rec(CoNode* x, int level=0){
          cout<<x<<"   "<<x->getnumber()<<" "<<level<<" ";
         if(x -> gettype() == Type::ZEROONE)cout<<"zeroone"<<endl;
        else cout<<"normal"<<endl;
        auto y = x -> get_head_of_list_of_children();
        int cnt = 0;
        while(y != nullptr){
            rec(y,level+1);
            y = y -> getnext();
        }
    }

    CoNode *Find_Lowest( CographRecognition::State & error){
        CoNode* y = new CoNode(Type::ZEROONE, 2);
        T -> add(y);
        CoNode *u, *w, *t;
        if(T -> getRoot() -> Marked_or_not() == Marked::UNMARKED){
            error = CographRecognition::State::COND_3;
            return y;
        }
        if(T -> getRoot() -> get_md() != T -> getRoot() -> get_d() - 1) {
            y = T->getRoot();
        }
            T -> getRoot() -> unmark();
            T -> getRoot() ->set_md(0);
            w = T -> getRoot();
            u = w;
        queue<CoNode*>q;
        Made_Queue_Of_Marked(T -> getRoot(), q);
        while(!q.empty()){
            u = q.front();
            q.pop();
            if(u -> Marked_or_not() != Marked::MARKED)continue;
            if(y -> getnumber() != 2){//1 or 2
                error = CographRecognition::State::COND_1;
                return y;
            }
            if(u -> getnumber() == 1){
                if(u->get_md() != u -> get_d() - 1) {
                    y = u;
                }
                if (u->getParent()->Marked_or_not() == Marked::MARKED) {//1 and 6
                    error = CographRecognition::State::COND_1;
                    return y;
                } else {
                    t = u->getParent()->getParent();
                }
            } else{
                y = u;
                t = u -> getParent();
            }
            u -> unmark();
            u ->set_md(0);
            while(t != w){
                if(t == T -> getRoot()){//4
                    error = CographRecognition::State::COND_4;
                    return y;
                }
                if(t -> Marked_or_not() != Marked::MARKED){//3 or 5 or 6
                    error = CographRecognition::State::COND_3;//!
                    return y;
                }
                if(t -> get_md() != t -> get_d() - 1){//2
                    error = CographRecognition::State::COND_2;
                    return y;
                }
                if(t -> getParent() -> Marked_or_not() == Marked::MARKED){//1
                    error = CographRecognition::State::COND_1;
                    return y;
                }
                t -> unmark();
                t ->set_md(0);
                t = t -> getParent() -> getParent();
            }
            w = u;
        }
        return w;
    }



    vector<CoNode *> get_were_marked(CoNode *u){
        auto x = u -> get_head_of_list_of_children();
        vector<CoNode *> a;
        while(x != nullptr && x->Marked_or_not() == Marked::MARKED_AND_UNMARKED){
            a.push_back(x);
            x = x -> getnext();
        }
        return a;
    }


    CoNode* get_last_from_children(CoNode* u){
        auto x = u -> get_head_of_list_of_children();
        while(x != nullptr && x->Marked_or_not() == Marked::MARKED_AND_UNMARKED){
            x = x -> getnext();
        }
        return x;
    }

    void Insert_x_to_CoTree(CoNode *u, CoNode *x){
        vector<CoNode*>a;
        int u_number = u -> getnumber();
        a = get_were_marked(u);
        if((a.size() == 1 && u_number == 0 )||(u -> get_d() - a.size() == 1 && u_number == 1)){
            CoNode* w = a[0];
            if(u_number == 1){
                w = get_last_from_children(u);
            }
            if(w -> gettype() == Type::VERTEX){
                CoNode* y = new CoNode(Type::ZEROONE, u_number ^ 1);
                T -> add(y);
                if(u_number == 0)u -> remove_were_marked();
                else u -> remove_were_not_marked();
                u->addchild(y);
                y -> addchild(x);
                y -> addchild(w);
            } else{
                w -> addchild(x);
            }
        } else{
            auto vec = u -> remove_were_marked();
            CoNode *y = new CoNode(Type::ZEROONE, u_number);
            T -> add(y);
            for(auto v : vec){
                y -> addchild(v);
            }
            if(u_number == 1){
                auto nxt = u -> getnext();
                auto prv = u -> getprev();
                if(prv != nullptr)prv ->setnext(y);
                if(nxt != nullptr)nxt -> setprev(y);
                y ->setprev(prv);
                y ->setnext(nxt);
                if(prv == nullptr && u -> getParent() != nullptr)u -> getParent() ->set_head_of_list_of_children(y);
                if(u -> getParent() != nullptr)y -> setParent(u -> getParent());
                else{
                    T ->setRoot(y);
                }
                CoNode* z = new CoNode(Type::ZEROONE, 0);
                T -> add(z);
                y -> addchild(z);
                z -> addchild(x);
                z -> addchild(u);
            } else{
                CoNode *z = new CoNode(Type::ZEROONE, 1);
                T -> add(z);
                u ->addchild(z);
                z ->addchild(x);
                z->addchild(y);
            }
        }
    }


    CographRecognition::State CographRecognition::Cograph_Recognition(NetworKit::Graph &graph){
        CoNode *R = new CoNode(Type::ZEROONE, 1);
        CoTree Tp(R);
        T = &Tp;
        T -> add(R);
        G = graph;
        vector<NetworKit::node>vertex;
        vector<CoNode*>covertex;
        map<NetworKit::node, int> pos;
        int cnt = 0;
        for(auto i : G.nodeRange()){
            vertex.push_back(i);
            pos[i] = cnt++;
            CoNode *C = new CoNode(Type::VERTEX,i);
            T -> add(C);
            covertex.push_back(C);
        }
        for(auto i : G.nodeRange()){
            vector<CoNode*> vec;
            for(auto u : G.neighborRange(i)){
                vec.push_back(covertex[pos[u]]);
            }
            covertex[pos[i]] -> setoutEdges(vec);
        }

        if(cnt == 0){
            T ->clear();
            return State::COGRAPH;
        }
        if(cnt == 1){
            R ->addchild(covertex[0]);
            T -> clear();
            return State::COGRAPH;
        }
        if(G.hasEdge(vertex[0], vertex[1])){
            R ->addchild(covertex[0]);
            R ->addchild(covertex[1]);
        } else {
            CoNode *N = new CoNode(Type::ZEROONE, 0);
            T->add(N);
            R->addchild(N);
            N->addchild(covertex[0]);
            N->addchild(covertex[1]);
        }
        covertex[0] -> add_to_graph();
        covertex[1] -> add_to_graph();
        CoNode* root = R;

        for(int i = 2; i < cnt; i++){
            root = T -> getRoot();
            Reset_All_CoNodes(root);
            Mark(covertex[i]);
            if(root -> Marked_or_not() == Marked::MARKED_AND_UNMARKED){//all nodes of T were marked and unmarked <=> R is marked and unmarked
                root -> addchild(covertex[i]);
            } else if(mark_ever_count == 0){
                if(root -> get_d() == 1){
                    (root -> get_head_of_list_of_children()) -> addchild(covertex[i]);
                } else{
                    CoNode *R1 = new CoNode(Type::ZEROONE, 1);
                    CoNode *R2 = new CoNode(Type::ZEROONE, 0);
                    T -> add(R1);
                    T -> add(R2);
                    R1 -> addchild(R2);
                    R2 -> addchild(root);
                    R2 -> addchild(covertex[i]);
                    T ->setRoot(R1);
                    root = R1;
                }
            } else{
                CographRecognition::State error = State::COGRAPH;
                CoNode* u = Find_Lowest(error);
                if(error != State::COGRAPH){
                    T -> clear();
                    return error;
                }
                Insert_x_to_CoTree(u, covertex[i]);
            }
            root = T -> getRoot();
            covertex[i] -> add_to_graph();
        }
        T -> clear();
        return State::COGRAPH;
    }


}