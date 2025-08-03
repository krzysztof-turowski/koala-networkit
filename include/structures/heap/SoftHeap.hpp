/**
 * SoftHeap.hpp
 * 
 *  Created on: 07.07.2025
 *      Author: Hubert Bernacki (hubert.bernacki@student.uj.edu.pl)
 */

#pragma once

#include <list>
#include <memory>
#include <optional>
#include <concepts>
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>

template<typename T>
concept SoftHeapElement = requires(T a){
    { a.key } -> std::convertible_to<int>;
    { a.ckey } -> std::convertible_to<int>;
    { a.corrupted } -> std::convertible_to<bool>;
};

template<SoftHeapElement T>
class SoftHeap {
public:
    // [TODO] implement other constructors
    SoftHeap() = delete;
    SoftHeap(T e, float eps = 0.3);
    SoftHeap(const SoftHeap&);
    SoftHeap(SoftHeap&&);
    // ~SoftHeap();
    SoftHeap& operator=(const SoftHeap&);
    SoftHeap& operator=(SoftHeap&&);

    template<SoftHeapElement R>
    friend SoftHeap<R> insert(SoftHeap<R>&&, R);
    T extractMin();
    void assertValidState(int expectedElements, bool checkSize = true, bool checkMinSuf = true);
    template<SoftHeapElement R>
    friend SoftHeap<R> meld(SoftHeap<R>&&, SoftHeap<R>&&);


    struct TreeNode {
        std::shared_ptr<TreeNode> left;
        std::shared_ptr<TreeNode> right;
        int ckey = 0;
        int rank = 0;
        int size = 1;
        std::list<T> l;
        
        TreeNode() {}
        TreeNode(T val) {
            ckey = val.key;
            l = std::list{val};
        }

        T pickElement();
        void sift();
        bool leaf();
        int corruptedCount();
    };

    struct ListNode {
        std::shared_ptr<TreeNode> tree;
        std::shared_ptr<ListNode> next;
        std::weak_ptr<ListNode> prev;
        std::weak_ptr<ListNode> sufMin;
        int rank = 0;
        bool isGuard = false;
    };
    std::shared_ptr<TreeNode> combine(std::shared_ptr<TreeNode> x, std::shared_ptr<TreeNode> y) {
        // std::cout << "combine" << std::endl;
        auto z = std::make_shared<TreeNode>();
        z->left = std::move(x);
        z->right = std::move(y);
        z->rank = z->left->rank + 1;
        if (z->rank <= r) {
            z->size = 1;
        } else {
            z->size = (3 * z->left->size + 1) / 2;
        }
        z->sift();
        return z;
    }
    void updateSuffixMin(std::shared_ptr<ListNode> t);
    void mergeInto(SoftHeap&& p);
    void repeatedCombine(int k);
    void insertTree(std::shared_ptr<ListNode> l1, std::shared_ptr<ListNode> l2);
    void removeTree(std::shared_ptr<ListNode> listNode);
    void assertValidListNode(std::shared_ptr<ListNode>, bool checkMinSuf = true);
    int assertValidTreeNode(std::shared_ptr<TreeNode>);
    int corruptedCount();
    std::shared_ptr<ListNode> first() {
        // std::cout << "first" << std::endl;
        return guard->next;
    }

    float eps;
    float r;
    int rank = 0;
    int insertCount = 0;
    std::shared_ptr<ListNode> guard;
};

// template<SoftHeapElement T>
// SoftHeap<T>::~SoftHeap(){

// }

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(T val, float eps) : eps(eps) {
    // std::cout << "SoftHeap::SoftHeap(T, float)" << std::endl;
    r = ceil(log(1 / eps) / log(2.0)) + 5;

    guard = std::make_shared<ListNode>();
    guard->isGuard = true;
    guard->rank = -1;

    auto node = std::make_shared<ListNode>();
    node->next = guard;
    node->prev = guard;
    node->sufMin = node;
    node->tree = std::make_shared<TreeNode>(val);

    guard->next = node;
    guard->prev = node;
    guard->sufMin = node;
}



// template<SoftHeapElement T>
// SoftHeap<T>::SoftHeap() {
//     // std::cout << "SoftHeap::SoftHeap()" << std::endl;
//     throw std::runtime_error("DEFAULT CONSTRUCTOR");
// }

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(const SoftHeap& other)
    : eps(other.eps), r(other.r), rank(other.rank), guard(other.guard)
{
    // std::cout << "SoftHeap::SoftHeap(const SoftHeap&)" << std::endl;
    // TODO It actually doesn't copy...
    throw std::runtime_error("COPY CONSTRUCTOR");
}

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(SoftHeap<T>&& other)
    : eps(other.eps), r(other.r), rank(other.rank), guard(other.guard)
{
    // std::cout << "SoftHeap::SoftHeap(SoftHeap&&)" << std::endl;
    other.guard.reset();
}

template<SoftHeapElement T>
SoftHeap<T>& SoftHeap<T>::operator=(const SoftHeap&) {
    // std::cout << "SoftHeap::operator=(const SoftHeap&)" << std::endl;
    
    throw std::runtime_error("COPY ASSIGNMENT");
}

template<SoftHeapElement T>
SoftHeap<T>& SoftHeap<T>::operator=(SoftHeap<T>&& other) {
    // std::cout << "SoftHeap::operator=(SoftHeap&&)" << std::endl;
    eps = other.eps;
    r = other.r;
    rank = other.rank;
    guard = other.guard;

    other.guard.reset();
    return *this;
}

template<SoftHeapElement T>
SoftHeap<T> insert(SoftHeap<T>&& p, T val) {
    // std::cout << "insert" << std::endl;
    auto pInserts = p.insertCount;
    auto sh = meld(std::move(p), std::move(SoftHeap<T>(val)));
    sh.insertCount = pInserts + 1;
    return sh;
}

template<SoftHeapElement T>
bool SoftHeap<T>::TreeNode::leaf() {
    // std::cout << "TreeNode::leaf" << std::endl;
    return !left && !right; 
}

template<SoftHeapElement T>
T SoftHeap<T>::TreeNode::pickElement() {
    // std::cout << "TreeNode::pickElement" << std::endl;
    if(l.empty()){
        throw std::runtime_error("EMPTY LIST POP");
    }
    T ret = l.front();
    ret.ckey = ckey;
    if (ret.ckey != ret.key) {
        ret.corrupted = true;
    }
    l.pop_front(); 
    return ret;
}

template<SoftHeapElement T>
void SoftHeap<T>::TreeNode::sift() {
    // std::cout << "TreeNode::sift" << std::endl;
    while (l.size() < size && !leaf()) {
        if(!left || (right && left->ckey > right->ckey)) {
            std::swap(left, right);
        }
        l.splice(l.end(), left->l);
        ckey = left->ckey;
        left->l = std::list<T>{};
        if (left->leaf()) {
            left.reset();
        } else {
            left->sift();
        }
    }

    assert(leaf() || (size <= l.size() && l.size() <= 3 * size));
}

template<SoftHeapElement T>
void SoftHeap<T>::removeTree(std::shared_ptr<SoftHeap<T>::ListNode> t) {
    // std::cout << "removeTree" << std::endl;
    t->prev.lock()->next = t->next;
    t->next->prev = t->prev.lock();
}

template<SoftHeapElement T>
void SoftHeap<T>::updateSuffixMin(std::shared_ptr<SoftHeap<T>::ListNode> t) {
    // std::cout << "updateSuffixMin" << std::endl;
    // assertValidState(0, false, false);
    while (!t->isGuard) {
        // std::cout << "while loop" << std::endl;
        // std::cout<< "next is guard: "<< t->next->isGuard << std::endl;
        assert(t->tree);
        if (t->next->isGuard 
            || t->tree->ckey <= t->next->sufMin.lock()->tree->ckey) {
            // std::cout << "if" << std::endl;
            t->sufMin = t;
        } else {
            // std::cout << "else" << std::endl;
            t->sufMin = t->next->sufMin;
        }
        
        // std::cout << "after if" << std::endl;
        t = t->prev.lock();
        // std::cout << "after prev.lock()" << std::endl;
    }
    // std::cout<<"exited updateSuffixMin" << std::endl;
}

template<SoftHeapElement T>
void SoftHeap<T>::repeatedCombine(int k) {
    // std::cout << "repeatedCombine" << std::endl;
    auto t = first();
    while (!t->next->isGuard) {
        if (t->rank == t->next->rank) {
            if (t->next->next->isGuard || t->rank != t->next->next->rank) {
                t->tree = combine(std::move(t->tree), std::move(t->next->tree));
                t->rank = t->tree->rank;
                removeTree(t->next);
            } else if (t->rank > k) {
                break;
            }
        } else if (t->rank > k) {
            // This break breaks probably because combine operation 
            // creates a heap with a larger rank so we need to keep combining... 
            // I think this else if was supposed to be inside the previous if...
            //
            // break;
            // std :: cout << "NO BREAKS" << std::endl;
        }
        t = t->next;
    }

    if (t->rank > rank) {
        rank = t->rank;
    }
    // assertValidState(0, false, false);
    assert(!t->isGuard);
    updateSuffixMin(t);
}


template<SoftHeapElement T>
SoftHeap<T> meld(SoftHeap<T>&& p, SoftHeap<T>&& q) {
    // std::cout << "meld" << std::endl;
    if(p.rank > q.rank) std::swap(p, q);
    int prank = p.rank;
    q.mergeInto(std::move(p));
    q.repeatedCombine(prank);
    return std::move(q);
}



template<SoftHeapElement T>
void SoftHeap<T>::insertTree(
    std::shared_ptr<SoftHeap<T>::ListNode> l1,
    std::shared_ptr<SoftHeap<T>::ListNode> l2
) {
    // std::cout << "insertTree" << std::endl;
    l1->next = l2;
    l2->prev.lock()->next = l1;
    l1->prev = l2->prev.lock();
    l2->prev = l1;
    // Why don't update every ptr ???
    // l2->prev = l1;
}
    

template<SoftHeapElement T>
void SoftHeap<T>::mergeInto(SoftHeap<T>&& p) {
    // std::cout << "mergeInto" << std::endl;
    if (p.rank > rank) throw std::runtime_error("MERGING LARGER TO SMALLER");

    auto t1 = p.first();
    auto t2 = first();

    while (t1 && !t1->isGuard) {
        while (t1->rank > t2->rank) {
            t2 = t2->next;
        }

        auto t1new = t1->next;
        insertTree(t1, t2);
        t1 = t1new;
    }
}

template<SoftHeapElement T>
T SoftHeap<T>::extractMin() {
    // std::cout << "extractMin" << std::endl;
    if (first()->isGuard) throw std::runtime_error("EXTRACT FROM EMPTY HEAP");
    auto t = first()->sufMin.lock();
    auto x = t->tree;
    auto e = x->pickElement();

    if (x->l.size() <= x->size / 2) {
        if (!x->leaf()) {
            x->sift();
            updateSuffixMin(t);
        }else if (x->l.empty()) {
            removeTree(t);
        }
    }
    // I think this is neccessary...
    // The problem is that it affects the complexity, potentially makes it
    // something like O(log (size))
    // updateSuffixMin(guard->prev.lock());
    
    // We may only need to call it if prev.sufMin == t
    // Also we could call it starting at the previous node...
    if (t->prev.lock()->sufMin.lock() == t) {
        updateSuffixMin(t->prev.lock());
    }
    return e;
}

template<SoftHeapElement T>
void SoftHeap<T>::assertValidListNode(std::shared_ptr<SoftHeap<T>::ListNode> node, bool checkMinSuf) {
    assert(node);
    assert(node->next);
    assert(node->prev.lock());
    assert(!checkMinSuf || node->isGuard || node->sufMin.lock());
    assert(!checkMinSuf || node->isGuard || !node->sufMin.lock()->isGuard);
    assert(node->next->prev.lock() == node);
    assert(node->prev.lock()->next == node);
    assert(node->isGuard || node->tree);

    // std :: cout << "List Node" << std::endl; 
    // std :: cout << "Is guard " << (node->isGuard) << std::endl;
    // std :: cout << "THIS " << ((long long)node.get() & 0xFFFFll) << std::endl;
    // std :: cout << "NEXT " << ((long long)node->next.get() & 0xFFFFll) << std::endl;
    // std :: cout << "PREV " << ((long long)node->prev.lock().get() & 0xFFFFll) << std::endl;
    // if (checkMinSuf) {
    //     std :: cout << "SUF MIN " << ((long long)node->sufMin.lock().get() & 0xFFFFll) << std::endl;
    // }
}

template<SoftHeapElement T>
int SoftHeap<T>::assertValidTreeNode(std::shared_ptr<SoftHeap<T>::TreeNode> node) {
    assert(node);
    int c = node->l.size();
    // std :: cout << "TREE NODE" << std::endl;
    // std :: cout << "RANK " << node->rank << std::endl;
    // std :: cout << "SIZE " << node->size << std::endl;
    // std :: cout << "LIST LENGTH " << node->l.size() << std::endl;

    for (auto e : node->l) {
        assert(e.key <= node->ckey);
    }

    if (node->left) {
        c += assertValidTreeNode(node->left);
        assert(node->rank == node->left->rank + 1);
    }
    if (node->right) {
        c += assertValidTreeNode(node->right);
        assert(node->rank == node->right->rank + 1);
    }

    return c;
}

template<SoftHeapElement T>
int SoftHeap<T>::TreeNode::corruptedCount() {
    int c = 0;
    for (T e : l) {
        if (e.key < ckey) {
            c += 1;
        }
    }

    if (left) {
        c += left->corruptedCount();
    }
    if (right) {
        c += right->corruptedCount();
    }
    return c;
}

template<SoftHeapElement T>
int SoftHeap<T>::corruptedCount() {
    auto node = first();
    int corrupted = 0;
    while (!node->isGuard) {
        corrupted += node->tree->corruptedCount();
        node = node->next;
    }
    return corrupted;
}

template<SoftHeapElement T>
void SoftHeap<T>::assertValidState(int expectedElements, bool checkSize, bool checkMinSuf) {
    assert(guard);
    assert(guard->isGuard);
    assertValidListNode(guard);
    auto node = guard->next;
    int elements = 0;
    
    while (node != guard) {
        assert(!node->isGuard);
        assertValidListNode(node, checkMinSuf);
        elements += assertValidTreeNode(node->tree);

        node = node->next;
    }
    // std :: cout << "elements: " << elements << std::endl; 
    assert(!checkSize || elements == expectedElements);
    // std :: cout << "CORRUPTED: " << corruptedCount() << std::endl;
    
    assert(corruptedCount() <= (eps) * insertCount);
}
/**
 * TODO 
 * 0. Top (without pop)
 * 1. Delete operation, lazy / delete from the linked list - maybe need SoftHeapElement to be ptr and have ->iterator so that we have access to the list...
 * 2. need access to all of the corrupted elements (!but only once! or maybe not idk), as elements cannot (? or maybe can uncorrupt themselves) it's easiest to just have a 
 *    list<T> corrupted to which we add elements when they become corrupted (what about corrupted -> not -> corrupted -> ... cycles???)
 * *3. assert that sufMin's link correctly
 * 
 * Remarks:
 * 1. If delete is not lazy it's not obvious what to do and if elements can uncorrupt themselves by some funny list change.
 * 2. Generally it maybe so that an element can uncorrupt itselves when changing lists but don't think so as the lists get concatenated => max only rises
 * 3. CKEY HANDLING IN SIFT SEEMS BS - ITS THE OPPOSITE OF WHAT SHOULD HAPPEN....
 */