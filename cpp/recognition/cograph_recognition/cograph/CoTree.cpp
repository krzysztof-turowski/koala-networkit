//
// Created by milana on 30.03.24.
//
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

        CoNode *head_of_list_of_children;
        CoNode *next, *prev;//in list of children of its parent
        CoNode *parent;
        std::vector<CoNode *>out_edges;//neighbours of cur vertex in G

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
        std::vector<CoNode *> RemoveWereMarked(){
            auto u = head_of_list_of_children;
            std::vector<CoNode *>vec;
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
        std::vector<CoNode *>save;
    public:
        CoNode *root;
        explicit CoTree(CoNode *root):root(root){
            save.push_back(root);
        }
        void Add(CoNode *x){
            save.push_back(x);
        }
        void Clear(){
            for(auto u : save){
                delete u;
            }
        }
    };
}