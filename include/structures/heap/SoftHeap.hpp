
#include <list>
#include <memory>
#include <optional>
#include <concepts>
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>
#include <variant>

template<typename T>
concept SoftHeapElement = requires(T a){
    { a->key } -> std::convertible_to<int>;
    { a->ckey } -> std::convertible_to<int>;
    { a->corrupted } -> std::convertible_to<bool>;
    { a->removed } -> std::convertible_to<bool>;
};

template<SoftHeapElement T>
class SoftHeap {
 public:
    explicit SoftHeap(float eps = 0.3);
    explicit SoftHeap(T e, float eps = 0.3);
    SoftHeap(const SoftHeap&);
    SoftHeap(SoftHeap&&);

    SoftHeap& operator=(const SoftHeap&);
    SoftHeap& operator=(SoftHeap&&);

    template<SoftHeapElement R>
    friend SoftHeap<R> insert(SoftHeap<R>&&, R);
    T extractMin();
    T extractMinInternal();
    std::variant<T, bool> lookupMin();
    std::variant<T, bool> lookupMinInternal();
    std::vector<T> corruptedElements();
    bool empty() { return elements <= 0; }

    void assertValidState(int expectedElements, bool checkSize = true, bool checkMinSuf = true);
    template<SoftHeapElement R>
    friend SoftHeap<R> meld(SoftHeap<R>&&, SoftHeap<R>&&);


    struct TreeNode {
        std::shared_ptr<TreeNode> left;
        std::shared_ptr<TreeNode> right;
        int ckey = 0;
        int rank = 0;
        int size = 1;
        std::list<T> corruptedList;
        std::list<T> originalList;

        TreeNode() {}
        explicit TreeNode(T val) {
            ckey = val->key;
            originalList = std::list{val};
        }

        void pushCorruptedElementsOnVector(std::vector<T>& vec);
        T pickElement(bool lookup = false);
        void sift();
        bool leaf();
        int corruptedCount();
        int listsSize();
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
    bool mergeInto(SoftHeap&& p);
    void repeatedCombine(int k);
    void insertTree(std::shared_ptr<ListNode> l1, std::shared_ptr<ListNode> l2);
    void removeTree(std::shared_ptr<ListNode> listNode);
    void assertValidListNode(std::shared_ptr<ListNode>, bool checkMinSuf = true);
    int assertValidTreeNode(std::shared_ptr<TreeNode>);
    int corruptedCount();
    std::shared_ptr<ListNode> first() {
        return guard->next;
    }
    std::list<T> allElementsAdded;

    int elements = 0;
    float eps;
    float r;
    int rank = 0;
    int insertCount = 0;
    std::shared_ptr<ListNode> guard;
};

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(float eps) : eps(eps), elements(0) {
    r = ceil(log(1 / eps) / log(2.0)) + 5;

    guard = std::make_shared<ListNode>();
    guard->isGuard = true;
    guard->rank = -1;

    guard->next = guard;
    guard->prev = guard;
    guard->sufMin = guard;
}

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(T val, float eps) : eps(eps), elements(1) {
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

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(const SoftHeap& other)
    : eps(other.eps), r(other.r), rank(other.rank), guard(other.guard), elements(elements)
{}

template<SoftHeapElement T>
SoftHeap<T>::SoftHeap(SoftHeap<T>&& other)
    : eps(other.eps), r(other.r), rank(other.rank), guard(other.guard), elements(elements) {
    other.guard.reset();
}

template<SoftHeapElement T>
SoftHeap<T>& SoftHeap<T>::operator=(const SoftHeap&) {
    throw std::runtime_error("COPY ASSIGNMENT");
}

template<SoftHeapElement T>
SoftHeap<T>& SoftHeap<T>::operator=(SoftHeap<T>&& other) {
    eps = other.eps;
    r = other.r;
    rank = other.rank;
    guard = other.guard;
    elements = other.elements;

    other.guard.reset();
    return *this;
}

template<SoftHeapElement T>
SoftHeap<T> insert(SoftHeap<T>&& p, T val) {
    auto pInserts = p.insertCount;
    auto sh = meld(std::move(p), std::move(SoftHeap<T>(val)));
    sh.insertCount = pInserts + 1;
    return sh;
}

template<SoftHeapElement T>
bool SoftHeap<T>::TreeNode::leaf() {
    return !left && !right;
}

template<SoftHeapElement T>
void SoftHeap<T>::TreeNode::pushCorruptedElementsOnVector(std::vector<T>& vec) {
    for ( T t : corruptedList ) {
        vec.push_back(t);
    }
    if (left) {
        left->pushCorruptedElementsOnVector(vec);
    }
    if (right) {
        right->pushCorruptedElementsOnVector(vec);
    }
}

template<SoftHeapElement T>
T SoftHeap<T>::TreeNode::pickElement(bool lookup) {
    if (originalList.empty() && corruptedList.empty()) {
        throw std::runtime_error("EMPTY LIST POP");
    }
    T ret;
    if (!originalList.empty()) {
        ret = originalList.front();
        if (!lookup) {
            originalList.pop_front();
        }
    } else {
        ret = corruptedList.front();
        if (!lookup) {
            corruptedList.pop_front();
        }
    }
    ret->ckey = ckey;
    if (ret->ckey != ret->key) {
        ret->corrupted = true;
    }
    return ret;
}

template<SoftHeapElement T>
int SoftHeap<T>::TreeNode::listsSize() {
    return originalList.size() + corruptedList.size();
}

template<SoftHeapElement T>
void SoftHeap<T>::TreeNode::sift() {
    while (listsSize() < size && !leaf()) {
        if (!left || (right && left->ckey > right->ckey)) {
            std::swap(left, right);
        }
        if (left->ckey != ckey) {
            for (T t : originalList) {
                t->corrupted = true;
                corruptedList.push_back(t);
            }
            originalList.clear();
        }
        originalList.splice(originalList.end(), left->originalList);
        corruptedList.splice(corruptedList.end(), left->corruptedList);
        ckey = left->ckey;
        left->originalList = std::list<T>{};
        left->corruptedList = std::list<T>{};

        if (left->leaf()) {
            left.reset();
        } else {
            left->sift();
        }
    }

    assert(leaf() || (size <= listsSize() && listsSize() <= 3 * size));
}

template<SoftHeapElement T>
void SoftHeap<T>::removeTree(std::shared_ptr<SoftHeap<T>::ListNode> t) {
    t->prev.lock()->next = t->next;
    t->next->prev = t->prev.lock();
}

template<SoftHeapElement T>
void SoftHeap<T>::updateSuffixMin(std::shared_ptr<SoftHeap<T>::ListNode> t) {
    while (!t->isGuard) {
        assert(t->tree);
        if (t->next->isGuard
            || t->tree->ckey <= t->next->sufMin.lock()->tree->ckey) {
            t->sufMin = t;
        } else {
            t->sufMin = t->next->sufMin;
        }
        t = t->prev.lock();
    }
}

template<SoftHeapElement T>
void SoftHeap<T>::repeatedCombine(int k) {
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
            //
            // break;
        }
        t = t->next;
    }

    if (t->rank > rank) {
        rank = t->rank;
    }
    assert(!t->isGuard);
    updateSuffixMin(t);
}

template<SoftHeapElement T>
SoftHeap<T> meld(SoftHeap<T>&& p, SoftHeap<T>&& q) {
    int newElements = p.elements + q.elements;
    if (std::holds_alternative<bool>(p.lookupMin())) return q;
    if (std::holds_alternative<bool>(q.lookupMin())) return p;
    if (p.rank > q.rank) std::swap(p, q);
    int prank = p.rank;
    q.mergeInto(std::move(p));
    q.repeatedCombine(prank);
    q.elements = newElements;
    return std::move(q);
}

template<SoftHeapElement T>
void SoftHeap<T>::insertTree(
    std::shared_ptr<SoftHeap<T>::ListNode> l1,
    std::shared_ptr<SoftHeap<T>::ListNode> l2
) {
    l1->next = l2;
    l2->prev.lock()->next = l1;
    l1->prev = l2->prev.lock();
    l2->prev = l1;
}

template<SoftHeapElement T>
bool SoftHeap<T>::mergeInto(SoftHeap<T>&& p) {
    if (p.rank > rank) throw std::runtime_error("MERGING LARGER TO SMALLER");

    auto t1 = p.first();
    auto t2 = first();

    while (t1 && !t1->isGuard) {
        while (t1->rank > t2->rank) {
            if (t2->isGuard) {
                return false;
            }
            t2 = t2->next;
        }

        auto t1new = t1->next;
        insertTree(t1, t2);
        t1 = t1new;
    }
    return true;
}

template<SoftHeapElement T>
std::vector<T> SoftHeap<T>::corruptedElements() {
    auto lnode = this->first();

    std::vector<T> ret;
    while (!lnode->isGuard) {
        lnode->tree->pushCorruptedElementsOnVector(ret);
        lnode = lnode->next;
    }
    return ret;
}

template<SoftHeapElement T>
std::variant<T, bool> SoftHeap<T>::lookupMinInternal() {
    if (first()->isGuard) return false;
    auto t = first()->sufMin.lock();
    auto x = t->tree;
    auto e = x->pickElement(true);
    return e;
}

template<SoftHeapElement T>
std::variant<T, bool> SoftHeap<T>::lookupMin() {
    while (true) {
        auto retVal = lookupMinInternal();
        if (std::holds_alternative<bool>(retVal)) return retVal;
        auto ret = std::get<T>(retVal);
        if (!ret->removed) {
            return ret;
        }
        extractMinInternal();
    }
}

template<SoftHeapElement T>
T SoftHeap<T>::extractMin() {
    while (true) {
        T ret = extractMinInternal();
        if (!ret->removed) {
            return ret;
        }
    }
}

template<SoftHeapElement T>
T SoftHeap<T>::extractMinInternal() {
    elements -= 1;
    assert(!first()->isGuard);
    if (first()->isGuard) throw std::runtime_error("EXTRACT FROM EMPTY HEAP");
    auto t = first()->sufMin.lock();
    auto x = t->tree;
    auto e = x->pickElement();

    if (x->listsSize() <= x->size / 2) {
        if (!x->leaf()) {
            x->sift();
            updateSuffixMin(t);
        } else if (x->listsSize() == 0) {
            removeTree(t);
        }
    }

    if (t->prev.lock()->sufMin.lock() == t) {
        updateSuffixMin(t->prev.lock());
    }
    return e;
}

template<SoftHeapElement T>
void SoftHeap<T>::assertValidListNode(std::shared_ptr<SoftHeap<T>::ListNode> node,
    bool checkMinSuf) {
    assert(node);
    assert(node->next);
    assert(node->prev.lock());
    assert(!checkMinSuf || node->isGuard || node->sufMin.lock());
    assert(!checkMinSuf || node->isGuard || !node->sufMin.lock()->isGuard);
    assert(node->next->prev.lock() == node);
    assert(node->prev.lock()->next == node);
    assert(node->isGuard || node->tree);
}

template<SoftHeapElement T>
int SoftHeap<T>::assertValidTreeNode(std::shared_ptr<SoftHeap<T>::TreeNode> node) {
    assert(node);
    int c = node->listsSize();

    for (auto e : node->originalList) {
        assert(e->key <= node->ckey);
    }
    for (auto e : node->corruptedList) {
        assert(e->key <= node->ckey);
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
    int corrupted_count = 0;
    int corrupted_count_original = std::count_if(originalList.begin(), originalList.end(),
        [this](T e) { return e->key < ckey; });
    int corrupted_count_corrupted = std::count_if(corruptedList.begin(), corruptedList.end(),
        [this](T e) { return e->key < ckey; });

    assert(corrupted_count_original == 0);
    corrupted_count = corrupted_count_corrupted + corrupted_count_corrupted;

    if (left) {
        corrupted_count += left->corruptedCount();
    }
    if (right) {
        corrupted_count += right->corruptedCount();
    }
    return corrupted_count;
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
}
