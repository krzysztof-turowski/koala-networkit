#pragma once

#include <set>
#include <list>
#include <vector>
#include <functional>
#include <iostream>

namespace Koala {

#define QUEUE_DEBUG 0

// TODO split these structures to different files / rename the file

/**
 * Implementation of a binary heap with ability to remove elements. Also allows for iterating
 * over all elements (not in order).
 * 
 * @tparam V type of stored element
*/
template<typename V>
class Heap {
public:
    using handle_type = size_t;
    using const_iterator = typename std::vector<V>::const_iterator;
    using iterator = const_iterator;
    

    Heap() {}

    bool empty() const { return heap.empty(); }

    void clear() { 
        heap.clear();
        handle.clear(); 
        heap_index.clear();
    }

    size_t size() const { return heap.size(); }

    handle_type push(V value) {
        handle_type ref = heap_index.size();
        size_t index = heap.size();
        heap.push_back(value);
        handle.push_back(ref);
        heap_index.push_back(index);
        push_up(index);

        return ref;
    }

    V top() const { return heap.front(); } 

    void pop() { erase(handle[0]); }

    void erase(handle_type ref) {
        size_t index = heap_index[ref];
        size_t last_index = heap.size() - 1;
        swap(index, last_index);
        heap.pop_back();
        handle.pop_back();
        if (index < heap.size()) {
            push_down(index);
            push_up(index);
        }
    }
    
    void for_elements(const std::function<void(V)>& handle) const {
        for (auto v : heap) 
            handle(v);
    }

    const_iterator begin() const { return heap.begin(); }

    const_iterator end() const { return heap.end(); }

private:
    std::vector<V> heap;
    std::vector<handle_type> handle;
    std::vector<size_t> heap_index;

    size_t parent_index(size_t index) { return (index - 1) / 2; }
    size_t left_son_index(size_t index) { return 2 * index + 1; }
    size_t right_son_index(size_t index) { return 2 * index + 2; }

    void push_up(size_t index) {
        while (index > 0) {
            size_t parent = parent_index(index);
            if (heap[index] < heap[parent]) {
                swap(index, parent);
                index = parent;
            } else {
                break;
            }
        }
    }

    void push_down(size_t index) {
        size_t minimum = index;
        size_t left_son = left_son_index(index);
        size_t right_son = right_son_index(index);;

        if (left_son < heap.size() && heap[left_son] < heap[minimum])
            minimum = left_son;
        if (right_son < heap.size() && heap[right_son] < heap[minimum])
            minimum = right_son;

        if (minimum != index) {
            swap(index, minimum);
            push_down(minimum);
        }
    }

    void swap(size_t index_1, size_t index_2) {
        std::swap(heap[index_1], heap[index_2]);
        std::swap(heap_index[handle[index_1]], heap_index[handle[index_2]]);
        std::swap(handle[index_1], handle[index_2]);
    }
};

/**
 * Implementation of pq1 from the paper 'An O(EV log V) algorithm for finding a maximal 
 * weighted mathcing in general graphs' by Galil, Micali, Gabow. 
 * Used to store elements and their priorities and additionally allows for 
 * decreasing all priorities by a given amount efficiently.
 * 
 * @tparam E type of elements
 * @tparam P type of priority values
*/
template<typename E, typename P>
class PriorityQueue1 {
public:
    /**
     * Creates a queue that stores elements between 0 and size - 1.
     * 
     * @param size the maximum value of elements that can be stored
    */
    PriorityQueue1(E size): 
        heap(),
        element_handle(size),
        in_heap(size, false),
        modified_priority(size, 0) {}

    /**
     * Insert an element with given priority
     * 
     * @param value the value of inserted element
     * @param priority the priority of inserted element
    */
    void insert(E value, P priority) {
        if (!in_heap[value]) {
            modified_priority[value] = priority + Delta;
            element_handle[value] = heap.push({ priority + Delta, value });
            in_heap[value] = true;
        } else if (priority < current_priority(value)) {
            remove(value);
            modified_priority[value] = priority + Delta;
            element_handle[value] = heap.push({ priority + Delta, value });
        }
    }

    /**
     * Remove the element with provided value
     * 
     * @param value the value of element to remove
    */
    void remove(E value) { 
        heap.erase(element_handle[value]); 
        in_heap[value] = false;
    }

    /**
     * Remove the minimum element
    */
    void remove_min() { 
        in_heap[heap.top().second] = false;
        heap.pop(); 
    }

    /**
     * Check if queue is empty
     * 
     * @return 'true' if queue is empty, 'false' otherwise
    */
    bool empty() const { return heap.empty(); }

    /**
     * Clear the queue
    */
    void clear() {
        for (auto [p, e] : heap) in_heap[e] = false;
        heap.clear();
        Delta = 0;
    }

    /**
     * Decrease all priorities of current elements by a provided value
     * 
     * @param delta the amount by which to decrease priorities
    */
    void decrease_all_priorities(P delta) { Delta += delta; }

    /**
     * Check the current priority of element
     * 
     * @param value the value of element to check
     * @returns current priority of the provided element
    */
    P current_priority(E value) const {
        return modified_priority[value] - Delta;
    }

    /**
     * Find the current minimum element
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P> find_min() const {
        auto min = heap.top();
        return { min.second, min.first - Delta };
    }

    /**
     * Call function for each element of queue. The function receives the value and priority of each element
     * 
     * @param handle the function to be called
    */
    void for_elements(const std::function<void(E,P)>& handle) {
        for (auto entry : heap) {
            handle(entry.second, entry.first - Delta);
        }
    }

// private:
public:
    P Delta = 0;
    using pair_heap = Heap<std::pair<P, E>>;
    pair_heap heap;
    std::vector<typename pair_heap::handle_type> element_handle;
    std::vector<bool> in_heap;
    std::vector<P> modified_priority;
};

/**
 * Implementation of concatenable queues. Stores elements with their priorities. 
 * Preserves the order of elements that is independent from their priorities.
 * Allows appending and inserting elements as well as concatenating queues and spliting them.
 * 
 * @tparam H type of head value stored in each queue (useful for identifying queues with find)
 * @tparam E type of element values
 * @tparam P type of priority values
*/
template<class H, class E, class P>
class ConcatenableQueue {
public:
    struct Node;
    using handle_type = Node*;

    H head;
    
    ConcatenableQueue(H head): head(head), root(nullptr) {}

    ConcatenableQueue(ConcatenableQueue&& other) {
        root = other.root;
        head = other.head;
        other.root = nullptr;
        fix_root();
        
        #if QUEUE_DEBUG
        check_consistency();
        #endif
    }

    ConcatenableQueue& operator=(ConcatenableQueue&& other) {
        if (this != &other) {
            root = other.root;
            head = other.head;
            other.root = nullptr;
            fix_root();
        }
        return *this;
    }

    ~ConcatenableQueue() {
        if (root != nullptr) {
            root->delete_children();
            delete root;
        }
    }

    // Delete copy constructor as it's not used
    // TODO implement copying maybe
    ConcatenableQueue(const ConcatenableQueue<H, E, P>& other) = delete;
    ConcatenableQueue& operator=(const ConcatenableQueue<H, E, P>& other) = delete;

    void swap(ConcatenableQueue<H, E, P>& other) {
        std::swap(root, other.root);
        std::swap(head, other.head);
        fix_root();
        other.fix_root();
    }

    /**
     * Append element to the end of queue
     * 
     * @param element value of appended element
     * @param priority priority of appended element
     * @return reference to inserted element, used for spliting and finding queues. 
     *         Preserved by splits and concatenations.
    */
    handle_type append(E element, P priority) {
        #if QUEUE_DEBUG
        std::cerr << this << " append " << element << std::endl;
        #endif

        if (empty()) {
            Node* new_node = new Node(this, element, priority);
            root = Node::from_child(this, new_node);
            return new_node;
        }

        handle_type last = last_element();
        return insert_after(last, element, priority);
    }

    /**
     * Insert a new element after provided one
     * 
     * @param ref reference to element after which the new one is to be inserted
     * @param element value of inserted element
     * @param priority priority of inserted element
     * @return reference to inserted element, used for spliting and finding queues. 
     *         Preserved by splits and concatenations.
    */
    handle_type insert_after(handle_type ref, E element, P priority) {
        Node* new_node = new Node(this, element, priority);

        add_child_after(ref->parent, ref, new_node);

        new_node->parent->update_min_until_root();

        return new_node;
    }

    /**
     * Insert a new element before provided one
     * 
     * @param ref reference to element before which the new one is to be inserted
     * @param element value of inserted element
     * @param priority priority of inserted element
     * @return reference to inserted element, used for spliting and finding queues. 
     *         Preserved by splits and concatenations.
    */
    handle_type insert_before(handle_type ref, E element, P priority) {
        Node* new_node = new Node(this, element, priority);

        add_child_before(ref->parent, ref, new_node);

        new_node->parent->update_min_until_root();

        return new_node;
    }

    /**
     * Check if queue is empty
     * 
     * @return 'true' if queue is empty, 'false' otherwise
    */
    bool empty() const { return root == nullptr; }

    /**
     * Remove an element
     * 
     * @param element_ref reference to the element to remove
    */
    void remove(handle_type element_ref) {     
        if (element_ref->parent == root && root->children_count() == 1) {
            delete root;
            delete element_ref;
            root = nullptr;
            return;
        }

        remove_node(element_ref);
    }

    /**
     * Find the current minimum element
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P> find_min() const {
        return { root->min_element, root->min_priority };
    }

    // TODO different memory management, maybe unique pointers?
    static ConcatenableQueue<H, E, P>* concat(
            ConcatenableQueue<H, E, P>&& left, ConcatenableQueue<H, E, P>&& right, 
            H head, bool update_min = true) {
        auto new_left = new ConcatenableQueue<H, E, P>(std::move(left));
        new_left->concat(std::move(right), head, update_min);
        return new_left;
    }

    void concat(ConcatenableQueue<H, E, P>&& other, H new_head, bool update_min = true) {
        if (other.empty()) { head = new_head; return; }
        if (empty()) {
            *this = std::move(other);
            head = new_head;
            return;
        }
        head = new_head;
        Node* left = root; cut(left);
        Node* right = other.root; cut(right);
        if (left->height == right->height) {
            root = Node::from_children(this, left, right);
        } else if (left->height > right->height) {
            while (left->height > right->height) left = left->rightmost_child();
            add_child_after(left->parent, left, right);
            if (update_min) left->parent->update_min_until_root();
        } else {
            while (right->height > left->height) right = right->leftmost_child();
            other.add_child_before(right->parent, right, left);
            if (update_min) right->parent->update_min_until_root();
            root = other.root;
        }
        fix_root();
        other.root = nullptr;
    }

    std::pair<ConcatenableQueue<H, E, P>*, ConcatenableQueue<H, E, P>*> 
    split(handle_type split_element, H head_left, H head_right) {
        Node* prev = split_element;
        Node* iter = split_element->parent;
        ConcatenableQueue<H, E, P>* left = new ConcatenableQueue<H, E, P>(head_left, split_element);
        ConcatenableQueue<H, E, P>* right = new ConcatenableQueue<H, E, P>(head_right);

        while (iter != nullptr) {
            int child_count = iter->children_count();
            int path = iter->child_index(prev);

            for (int i = path - 1; i >= 0; -- i) { 
                left = concat(iter->children[i], std::move(*left), head_left, false);
            }

            for (int i = path + 1; i < child_count; ++ i) {
                right = concat(std::move(*right), iter->children[i], head_right, false);
            }

            if (prev != split_element) delete prev;
            prev = iter;
            iter = iter->parent;
        }

        if (prev != split_element) delete prev;

        root = nullptr;
        left->last_element()->parent->update_min_until_root();
        if (!right->empty()) right->first_element()->parent->update_min_until_root();

        return { left, right };
    }

    void for_each(const std::function<void(E, P)>& handle) {
        if (root != nullptr) root->for_each(handle);
    }

    void print(std::ostream& out = std::cerr) {
        out << head << " | ";
        if (root != nullptr) {
            out << find_min().first << " : ";
            root->print(out);
        }
        out << std::endl;
    }

    void print_elements(std::ostream& out = std::cerr) {
        if (empty())
            out << head << " : ";
        else 
            out << head << " | " << find_min().first << " : ";
        out << "{";
        std::vector<E> elements;
        for_each([&elements] (E element, P p) { elements.push_back(element); });
        for (int i = 0; i < elements.size(); ++ i) {
            out << elements[i];
            if (i != elements.size() - 1) out << ", ";
        }
        out << "}\n";
    }

    void check_consistency() {
        if (root == nullptr) return;
        if (root->parent != nullptr) {
            std::cerr << "root's parent is not null\n";
            exit(1);
        }
        if (root->queue != this) {
            std::cerr << "root points at wrong queue\n";
            exit(1);
        }
        root->check_children();
    }

    struct Node {

        ConcatenableQueue* find_queue() {
            Node* iter = this;
            while (iter->parent != nullptr) iter = iter->parent;
            return iter->queue;
        }
    
    private:

        ConcatenableQueue<H, E, P>* queue; // Correct in the root
        Node* parent;
        Node* children[3];
        int height;
        E min_element;
        P min_priority;

        Node(ConcatenableQueue<H, E, P>* queue, E element, P priority): 
            queue(queue),
            parent(nullptr),
            children{nullptr, nullptr, nullptr},
            height(0),
            min_element(element),
            min_priority(priority) {}

        static Node* from_child(ConcatenableQueue<H, E, P>* queue, Node* child) {
            Node* node = new Node(queue, child->min_element, child->min_priority);
            node->children[0] = child;
            node->height = child->height + 1;
            child->parent = node;
            return node;
        }

        static Node* from_children(ConcatenableQueue<H, E, P>* queue, Node* left, Node* right) {
            Node* node = new Node(queue, left->min_element, left->min_priority);
            node->children[0] = left;
            node->children[1] = right;
            if (right->min_priority < node->min_priority) {
                node->min_priority = right->min_priority;
                node->min_element = right->min_element;
            }
            node->height = left->height + 1;
            left->parent = node;
            right->parent = node;
            return node;
        }

        int children_count() {
            if (children[0] == nullptr) return 0;
            if (children[1] == nullptr) return 1;
            if (children[2] == nullptr) return 2;
            return 3;
        }

        void update_min() {
            min_element = children[0]->min_element;
            min_priority = children[0]->min_priority;
            for (int i = 1; i < 3; ++ i) {
                if (children[i] == nullptr) break;
                if (children[i]->min_priority < min_priority) {
                    min_priority = children[i]->min_priority;
                    min_element = children[i]->min_element;
                }
            }
        }

        void update_min_until_root() {
            update_min();
            
            if (parent != nullptr) { 
                parent->update_min_until_root();
            } 
        }
    
        int child_index(Node* child) {
            for (int i = 0; i < 3; ++ i) 
                if (children[i] == child)
                    return i;
            return -1;
        }

        Node* rightmost_child() {
            return children[children_count()-1];
        }

        Node* leftmost_child() {
            return children[0];
        }

        Node* left_sibling(Node* child) {
            int index = child_index(child);
            return index == 0 ? nullptr : children[index - 1];
        }

        Node* right_sibling(Node* child) {
            int child_count = children_count();
            int index = child_index(child);
            return index == child_count - 1 ? nullptr : children[index + 1];
        }

        Node* sibling(Node* child) {
            Node* left = left_sibling(child);
            return left == nullptr ? right_sibling(child) : left;
        }

        bool is_leaf() {
            return children[0] == nullptr;
        }

        void check_children() {
            if (parent == this) {
                std::cerr << "parent link points to the node itself\n";
                exit(1);
            }
            for (int i = 0; i < 3; ++ i) {
                if (children[i]) {
                    if (children[i]->parent != this) {
                        std::cerr << "parent link wrong "; print(); std::cerr << std::endl;
                        exit(1);
                    }
                    children[i]->check_children();
                }
            }
        }

        void for_each(const std::function<void(E, P)>& handle) {
            if (is_leaf()) {
                handle(min_element, min_priority);
                return;
            }
            for (int i = 0; i < 3; ++ i) {
                if (children[i] == nullptr) return;
                children[i]->for_each(handle);
            }
        }

        void print(std::ostream& out = std::cerr) {
            out << "(";
            if (is_leaf()) {
                out << min_element;
            } else {
                out << children_count() << ":";
                for (int i = 0; i < 3; ++ i) {
                    if (children[i] == nullptr) break;
                    children[i]->print(out);
                }
            }
            out << ")";
        }

        void delete_children() {
            for (int i = 0; i < 3; ++ i)
                if (children[i] != nullptr) {
                    children[i]->delete_children();
                    delete children[i];
                }
        }

        friend class ConcatenableQueue<H, E, P>;
    };

private:
    Node* root;

    ConcatenableQueue(H head, E element, P priority): 
        head(head), 
        root(new Node(element, priority)) {}

    ConcatenableQueue(H head, Node* element): 
            head(head) {
        if (element->height == 0) {
            root = Node::from_child(this, element);
            element->parent = root;
        } else {
            root = element;
            root->parent = nullptr;
        }

        fix_root();

        #if QUEUE_DEBUG
        check_consistency();
        #endif
    }

    static ConcatenableQueue<H, E, P>* concat(
            ConcatenableQueue<H, E, P>&& left, ConcatenableQueue<H, E, P>::Node* right, 
            H head, bool update_min = true) {
        ConcatenableQueue<H, E, P> right_que(head, right);
        return concat(std::move(left), std::move(right_que), head, update_min);
    }

    static ConcatenableQueue<H, E, P>* concat(
                ConcatenableQueue<H, E, P>::Node* left, ConcatenableQueue<H, E, P>&& right, 
            H head, bool update_min = true) {
        ConcatenableQueue<H, E, P> left_que(head, left);
        return concat(std::move(left_que), std::move(right), head, update_min);
    }

    void fix_root() { 
        if (root != nullptr) { 
            root->parent = nullptr; 
            root->queue = this; 
        } 
    }

    void insert_child(Node** children, int index, int child_count, Node* child) {
        for (int i = child_count; i > index; -- i)
            children[i] = children[i-1];
        children[index] = child;
    }

    void make_new_root(Node* a, Node* b) {
        root = Node::from_children(this, a, b);
    }

    void add_child_after(Node* target, Node* after, Node* child) {
        if (target == nullptr) {
            make_new_root(after, child);
            return;
        }
        
        int after_index = target->child_index(after);
        insert_child_at_index(target, after_index + 1, child);
    }

    void add_child_before(Node* target, Node* before, Node* child) {
        if (target == nullptr) {
            make_new_root(child, before);
            return;
        }
        
        int before_index = target->child_index(before);
        insert_child_at_index(target, before_index, child);
    }

    void insert_child_at_index(Node* target, int index, Node* child) {
        int child_count = target->children_count();
        if (child_count < 3) {
            insert_child(target->children, index, child_count, child);
            child->parent = target;
        } else {
            Node* four_children[4] = { target->children[0], target->children[1], 
                                       target->children[2], nullptr };
            insert_child(four_children, index, child_count, child);
            target->children[0] = four_children[0]; four_children[0]->parent = target;
            target->children[1] = four_children[1]; four_children[1]->parent = target;
            target->children[2] = nullptr;
            Node* new_node = Node::from_children(this, four_children[2], four_children[3]);
            new_node->parent = target->parent;
            target->update_min();
            add_child_after(target->parent, target, new_node);
        }
    }

    void remove_child(Node* from, Node* child) {
        int children = from->children_count();
        int child_index = from->child_index(child);
        for (int i = child_index; i < children - 1; ++ i)
            from->children[i] = from->children[i+1];
        from->children[children - 1] = nullptr;
        from->update_min();
    }

    handle_type first_element() {
        Node* iter = root;
        while (!iter->is_leaf()) {
            iter = iter->leftmost_child();
        }
        return iter;
    }

    handle_type last_element() {
        Node* iter = root;
        while (!iter->is_leaf()) {
            iter = iter->rightmost_child();
        }
        return iter;
    }

    void cut(Node*& node) {
        if (node->children_count() > 1) return;
        node = node->children[0]; 
        delete node->parent;
        node->parent = nullptr;
    }

    void remove_node(Node* node) {
        if (node->parent->children_count() == 3) {
            remove_child(node->parent, node);
            node->parent->update_min_until_root();
            delete node;
            return;
        }
        if (node->parent == root) {
            remove_child(root, node);
            if (root->height > 1) {
                Node* old_root = root;
                root = root->children[0];
                delete old_root;
                fix_root();
            }
            delete node;
            return;
        }
        Node* parent = node->parent;
        Node* sibling = parent->sibling(node);
        Node* left_uncle = parent->parent->left_sibling(parent);
        if (left_uncle != nullptr) {
            if (left_uncle->children_count() == 2) {
                left_uncle->children[2] = sibling;
                sibling->parent = left_uncle;
                left_uncle->update_min();
                remove_node(parent);
            } else {
                Node* move = left_uncle->children[2];
                remove_child(left_uncle, move);
                insert_child(parent->children, 0, 2, move);
                move->parent = parent;
                remove_child(parent, node);
                parent->update_min_until_root();
            }
            delete node;
            return;
        }
        Node* right_uncle = parent->parent->right_sibling(parent);
        if (right_uncle->children_count() == 2) {
            insert_child(right_uncle->children, 0, 2, sibling);
            sibling->parent = right_uncle;
            right_uncle->update_min();
            remove_node(parent);
        } else {
            Node* move = right_uncle->children[0];
            remove_child(right_uncle, move);
            parent->children[2] = move;
            move->parent = parent;
            remove_child(parent, node);
            parent->update_min_until_root();
        }
        delete node;
        return;
    }
};

/**
 * Implementation of pq2 from the paper 'An O(EV log V) algorithm for finding a maximal 
 * weighted mathcing in general graphs' by Galil, Micali, Gabow. 
 * Used to store elements and their priorities.
 * The elements are divided into groups that can be either active or nonactive.
 * Allows for changing the status of a group, changing priorities off all elements in active groups,
 * splitting groups
 * 
 * @tparam E type of elements
 * @tparam V type of values
 * @tparam P type of priority values
*/
template<typename E, typename V, typename P>
class PriorityQueue2 {
public:
    class Group { 
    public:

        bool empty() const {
            return elements.empty();
        }

        std::pair<E, P> find_min() const {
            auto [value, priority] = elements.find_min();
            return {value, priority - Delta_group};
        }

        bool is_active() const {
            return active;
        }

    private:
        Group(bool active, P Delta_last, P Delta_group): 
            elements(this), active(active), Delta_last(Delta_last), Delta_group(Delta_group) {}

        Group(ConcatenableQueue<Group*, E, P>&& elements_, 
                bool active, P Delta_last, P Delta_group): 
            elements(std::move(elements_)),
            active(active), Delta_last(Delta_last), Delta_group(Delta_group) {

            elements.head = this;
        }

        void update_Delta(P Delta) {
            if (active) {
                Delta_group = Delta_group + Delta - Delta_last;
            }
            Delta_last = Delta;
        }

        ConcatenableQueue<Group*, E, P> elements;
        bool active;
        P Delta_last;
        P Delta_group;

        friend class PriorityQueue2;
    };
    
    /**
     * Creates a queue that stores elements between 0 and size - 1.
     * 
     * @param size the maximum value of elements that can be stored
     * @param no_element placeholder for lack of element
     * @param no_value placeholder value
     * @param infinite_priority upper limit for priority
    */
    PriorityQueue2(E size, E no_element, V no_value, P infinite_priority): 
        no_element(no_element),
        no_value(no_value),
        infinite_priority(infinite_priority),
        group_minima(size), 
        elements(size, nullptr) {}

    /**
     * Append an element with given priority to the end group
     * 
     * @param element the value of inserted element
     * @param priority the priority of inserted element
     * @param group the group to which the new element belongs
    */
    void append(E element, P priority, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group_minimum(group);

        elements[element] = group->elements.append(element, priority + group->Delta_group);

        // Check if the new element is the minimum in group
        if (group->active && group->find_min().first == element) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(element, priority);
        }
    }

    /**
     * Insert an element with given priority into group before element
     * 
     * @param element the value of inserted element
     * @param priority the priority of inserted element
     * @param before the value of element before which the new element is inserted
     * @param group the group to which the new element belongs
    */
    void insert_before(E element, P priority, E before, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group_minimum(group);

        elements[element] = group->elements.insert_before(elements[before], element, priority + group->Delta_group);
        
        // Check if the new element is the minimum in group
        if (group->active && group->find_min().first == element) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(element, priority);
        }
    }

    /**
     * Remove an element
     * 
     * @param element the value of the element to remove
    */
    void remove(E element) { 
        Group* group = elements[element]->find_queue()->head;
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group_minimum(group);
        
        group->elements.remove(elements[element]);

        // Check if the removed element was the minimum in group
        if (last_min == element && group->active) {
            group_minima.remove(element);
            if (!group->empty()) {
                auto [new_min, new_min_priority] = group->find_min();
                group_minima.insert(new_min, new_min_priority);
            }
        }
    }

    /**
     * Check if there are any elements in active groups
     * 
     * @return 'true' if there is an active nonempty group, 'false' otherwise
    */
    bool has_active_elements() const {
        return !group_minima.empty();
    }

    /**
     * Find the current minimum element in any active group
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P> find_min() const {
        return group_minima.find_min();
    }

    /**
     * Decrease all priorities of current elements in active groups by a provided value
     * 
     * @param delta the amount by which to decrease priorities
    */
    void decrease_all_priorities(P delta) { 
        Delta += delta; 
        group_minima.decrease_all_priorities(delta);
    }

    /**
     * Create a new group
     * 
     * @param active status of new group
     * @return reference to the new group
    */
    Group* new_group(bool active) {
        return new Group(active, Delta, 0);
    }


    /**
     * Delete a group
     * 
     * @param group the group to be deleted
    */
    void delete_group(Group* group) {
        if (group->active && !group->empty()) {
            group_minima.remove(group->find_min().first);
        }
        delete group;
    }

    /**
     * Change the status of a group
     * 
     * @param group group whose status is to be changed
     * @param active the new status of the group
    */
    void change_status(Group* group, bool active) {
        if (group->active == active) return;
        group->update_Delta(Delta);
        group->active = active;

        if (active && !group->empty()) {
            auto [v, p] = group->find_min();
            group_minima.insert(v, p);
        } else if (!active && !group->empty()) {
            auto [v, p] = group->find_min();
            group_minima.remove(v);
        }
    }

    /**
     * Split a given group acording to the given element
     * 
     * @param group the group to be split
     * @param element the value of the element according to which the groups are to be split
     * @return a pair containing the two new groups, the first contains elements from beginning of group
     *      until and containing the split element, the second one contains elements after the element
    */
    std::pair<Group*, Group*> split_group(Group* group, E element) {
        group->update_Delta(Delta);

        if (group->active && !group->empty()) {
            auto [v, p] = group->find_min();
            group_minima.remove(v);
        }

        auto [elements_left, elements_right] = group->elements.split(elements[element], group, group);

        Group* group_left  = new Group(std::move(*elements_left),
                                group->active, group->Delta_last, group->Delta_group);
        Group* group_right = new Group(std::move(*elements_right),
                                group->active, group->Delta_last, group->Delta_group);

        if (group->active) {
            if(!group->empty()) {
                auto [v, p] = group->elements.find_min();
                group_minima.remove(v);
            }
            if (!group_left->empty()) {
                auto [v, p] = group_left->elements.find_min();
                group_minima.insert(v, p);
            }
            if (!group_right->empty()) {
                auto [v, p] = group_right->elements.find_min();
                group_minima.insert(v, p);
            }
        }
        delete group;
        return {group_left, group_right};
    }

    void shift_group(Group* group, E element) {
        group->update_Delta(Delta);

        auto [elements_left, elements_right] = group->elements.split(elements[element], group, group);

        elements_right->concat(std::move(*elements_left), group);
        group->elements = std::move(*elements_right);
        delete elements_left;
        delete elements_right;
    }

    void for_each_in_group(Group* group, const std::function<void(E, P)>& handle) {
        group->update_Delta(Delta);
        group->elements.for_each([&handle, group] (E e, P p) {
            handle(e, p - group->Delta_group);
        });
    }

    void for_group_minima(const std::function<void(E, P)>& handle) {
        group_minima.for_elements(handle);
    }

private:

    std::pair<E, P> group_minimum(Group* group) {
        return group->empty() ? std::make_pair(no_element, infinite_priority) : group->find_min();
    }

    E no_element;
    V no_value;
    P infinite_priority;
    PriorityQueue1<E, P> group_minima;
    std::vector<typename ConcatenableQueue<Group*, E, P>::handle_type> elements;
    
    P Delta = 0;
};

/**
 * Implementation of the union-find data structure operating on elements from a provided 
 * interval. Allows for linking sets of two provided elements.
 * Each set is associated with a value that can be set when linking two sets and that is returned
 * when calling the find method.
 * 
 * @tparam I type of elmenents
 * @tparam R type of the value assigned to each set
*/
template <typename I, typename R>
class UnionFind {
public:
    
    /**
     * Creates a union-find data structure for storing elements between 0 and size - 1
     * 
     * @param size upper bound on the elements that can be stored
    */
    UnionFind(I size): size(size), root(size), parent(size), rank(size) {}

    /**
     * Resets the data structure for a provided element. After calling the method the element
     * is in a set associated with a provided value. This is set might contain more than the 
     * provided element. In order for the data structure to represent only singletons the method 
     * has to be called for all elements that have been previously linked.
     * 
     * This is done individually so that operations can be done on an arbitrary subset of elements
     * between 0 and size - 1 without having to reset all size values.
     * 
     * @param element value of the element to be reset
     * @param value value that is to be associated with the element's new set
    */
    void reset(I element, R value) {
        parent[element] = element;
        root[element] = value;
        rank[element] = 0;
    }

    /**
     * Links sets of two provided elements.
     * 
     * @param x first element
     * @param y second element
     * @param value the value that is to associated with the created set
    */
    void link(I x, I y, R value) {
        I rx = find_root(x);
        I ry = find_root(y);

        if (rank[rx] < rank[ry]) {
            parent[rx] = ry;
            root[ry] = value;
        } else {
            parent[ry] = rx;
            root[rx] = value;
            if (rank[rx] == rank[ry]) {
                rank[rx] ++;
            }
        } 
    }

    /**
     * Returns the value associated with the set that the provided element belongs to.
     * 
     * @param element value of the element
    */
    R find(I element) const {
        return root[find_root(element)];
    }

private:
    I size;
    std::vector<R> root;
    std::vector<I> parent;
    std::vector<int> rank;

    I find_root(I x) {
        return parent[x] == x ? x : (parent[x] = find_root(parent[x]));
    }
};

/**
 * Implementation of the split-findmin data structure based on the list splitting algorithm
 * as described in the paper 'A Scaling Algorithm for Weighted Matching on General Graphs'
 * by Harold N. Gabow.
 * 
 * The data structure deals with elements between 0 and and a provided size. Every element can be in
 * at most one list and has an associated cost. The cost of a list is the smallest cost of any of 
 * it's elements. Each list has an id.
 * The structure allows for decreasing a cost a specific element to a given value, splitting
 * lists on provided element and finding the costs of lists and the elements that achieve the
 * minimum cost.
*/
template<typename T, typename H, typename K, typename V>
class SplitFindMin {
public: 
    class List;
    /**
     * Creates a split-findmin data structure with level alpha(m, n)
     * 
     * @param size maximum value of a stored element
     * @param infinity default value for a key
     * @param empty_val default value for a value associated with a key
     * @param n parameter used in calculating the level 
     * @param m parameter used in calculating the level
    */
    SplitFindMin(T size, K infinity, V empty_val, int n, int m):
        SplitFindMin(size, infinity, empty_val, alpha(m, n)) {}

    /**
     * Creates a list from a provided lists of elements with a given id.
     * All elements have default key values.
     * 
     * @param nodes elements to create a list of
     * @param id id associated with the list
     * @returns pointer to the new list
    */
    List* init(const std::list<T>& nodes, H id) {
        List* L = new List;
        L->id = id;
        L->i = max_i;
        L->sublist = nullptr;
        L->min_key = infinity;
        L->min_val = empty_val;

        initialize_head(L, nodes, id, max_i);

        return L;
    }

    /**
     * Returns the id of the list containing the provided element
     * 
     * @param element value of the element whose list's id is to be returned
     * @returns the id of the list containing the provided element
    */
    H list(T element) {
        return find_list(element, max_i)->id;
    }

    /**
     * Splits the list containing the provided element into two list - one containing all the
     * elements up to the provided one and a second one containg the remaining elements. The two lists
     * are given a new id.
     * 
     * @param element the element on whose position to split it's list
     * @param id1 the new id of the list containg elements up to the provided one
     * @param id2 the new id of the list containg remaining elements
     * @returns pair of two lists resulting from the split
    */
    std::pair<List*, List*> split(T element, H id1, H id2) {
        return split(element, id1, id2, max_i);
    }

    /**
     * Updates the key of the element if the new one is lower
     * 
     * @param element the element whose key is to be updated
     * @param key value of the key
     * @param value value associated with the key
    */
    void decreaseKey(T element, K key, V value) {
        decrease_key(element, key, value, max_i);
    }

    /**
     * Returns the minimum key of an element of the list and it's associated value
     * 
     * @param L list whose minimum key is returned
     * @returns pair containing the minimum key in the list and it's associated value
    */
    std::pair<K, V> findMin(List* L) {
        return {L->min_key, L->min_val};
    }

    /**
     * Returns the current key of an element and it's associated value
     * 
     * @param element elements for which current key and it's value is returned
     * @returns pair containing the current key of the provided element and it's associated value
    */
    std::pair<K, V> currentKey(T element) {
        return {key[max_i][element], val[max_i][element]};
    }

    /**
     * Deletes the list
     * 
     * @param L the list to be deleted
    */
    void deleteList(List* L) {
        reset_keys(L);
        delete L;
    }

    struct Sublist {
        int level;
        bool head;
        List* list;
        List* elements;
        std::list<Sublist*>::iterator sublist_it;

        ~Sublist() { delete elements; }
    };

    class List {
    public:
        void clear() {
            head_singletons.clear();
            for (auto s : head) delete s;
            head.clear();
            for (auto s : tail) delete s;
            tail.clear();
            tail_singletons.clear();
            nodes.clear();
        }

        H id;
        int i;
        K min_key;
        V min_val;
        std::list<T> nodes;
        std::list<Sublist*> head, tail;
        std::list<T> head_singletons, tail_singletons;
        Sublist* sublist;

        ~List() {
            clear();
        }
    };

private:

    static constexpr int inf_size = 1000000000;
    static constexpr int A2[4] = {2, 4, 16, 65536};

    static int A(int i, int j) {
        if (j == 0) return 2;
        if (i == 1) return j < 30 ? (1 << j) : inf_size;
        if (i == 2) return j < 4 ? A2[j] : inf_size;
        if (i == 3) return j == 1 ? 16 : inf_size;
        return inf_size;
    }

    static int a(int i, int n) {
        int j = -1;
        while (2 * A(i, j + 1) <= n) j ++;
        return j;
    }

    static int alpha(int m, int n) {
        int i = 1;
        while (A(i, m / n) < n) i ++;
        return i;
    }

    SplitFindMin(T size, K infinity, V empty_val, int level):
        max_i(level),
        size(size), 
        infinity(infinity), 
        empty_val(empty_val),
        e(level + 1),
        key(level + 1),
        val(level + 1),
        element_list(level + 1),
        list_it(level + 1),
        superelement_nodes(level + 1) {
            for (int i = 0; i <= level; ++ i) {
                e[i] = std::vector<T>(size, -1);
                element_list[i] = std::vector<List*>(size, nullptr);
                key[i] = std::vector<K>(size, infinity),
                val[i] = std::vector<V>(size, empty_val),
                list_it[i] = std::vector<typename std::list<T>::iterator>(size);
                superelement_nodes[i] = std::vector<std::list<T>>(size);
            }
        }
    
    int max_i;
    T size;
    K infinity;
    V empty_val;
    std::vector<std::vector<T>> e;
    std::vector<std::vector<K>> key;
    std::vector<std::vector<V>> val;
    std::vector<std::vector<List*>> element_list;
    std::vector<std::vector<typename std::list<T>::iterator>> list_it;
    std::vector<std::vector<std::list<T>>> superelement_nodes;

    void initialize_head(List* L, const std::list<T>& nodes, H id, int i) {
        auto nodes_cpy = nodes;
        for (auto it = nodes_cpy.begin(); it != nodes_cpy.end(); ++ it)
            list_it[i][*it] = it;
        L->nodes.splice(L->nodes.begin(), std::move(nodes_cpy));

        int j = a(i, nodes.size()) + 1;
        int remaining = nodes.size();
        auto it = nodes.rbegin();
        Sublist* sublist = nullptr;
        std::list<T> sublist_elements;

        do {
            int prev_j = j;
            while (j >= 0 && 2 * A(i, j) > remaining) -- j;

            if (j != prev_j) {
                if (sublist != nullptr) {
                    sublist->elements = new List;
                    sublist->elements->sublist = sublist;
                    sublist->elements->i = i - 1;
                    sublist->elements->id = id;
                    sublist->elements->min_key = infinity;
                    sublist->elements->min_val = empty_val;

                    initialize_head(sublist->elements, sublist_elements, id, i - 1);
                    sublist_elements.clear();
                    
                    L->head.push_front(sublist);
                    sublist->sublist_it = L->head.begin();
                }
                if (j == -1) {
                    sublist = nullptr;
                } else {
                    sublist = new Sublist;
                    sublist->list = L;
                    sublist->level = j;
                    sublist->head = true;
                }
            }    

            if (j == -1) {
                T singleton = *it;
                e[i][singleton] = -1;
                element_list[i][singleton] = L;

                L->head_singletons.push_front(singleton);
                it ++;
            } else {
                int size = 2 * A(i, j);
                T superelement = *it;
                sublist_elements.push_front(superelement);
                key[i-1][superelement] = infinity;
                val[i-1][superelement] = empty_val;

                std::list<T> se_nodes;
                for (int j = 0; j < size; ++ j) {
                    se_nodes.push_front(*it);
                    it ++;
                }

                for (auto n : se_nodes) {
                    e[i][n] = superelement;
                    update_key_at_level(superelement, i - 1, key[i][n], val[i][n]);
                }

                superelement_nodes[i][superelement] = std::move(se_nodes);

                remaining -= size;
            }
        } while (it != nodes.rend());

        if (sublist != nullptr) {
            sublist->elements = new List;
            sublist->elements->sublist = sublist;
            sublist->elements->i = i - 1;
            sublist->elements->id = id;
            sublist->elements->min_key = infinity;
            sublist->elements->min_val = empty_val;

            initialize_head(sublist->elements, sublist_elements, id, i - 1);
            
            L->head.push_front(sublist);
            sublist->sublist_it = L->head.begin();
        }

        for (auto n : nodes) update_key(L, key[i][n], val[i][n]);
    }

    void initialize_tail(List* L, const std::list<T>& nodes, H id, int i) {
        auto nodes_cpy = nodes;
        for (auto it = nodes_cpy.begin(); it != nodes_cpy.end(); ++ it)
            list_it[i][*it] = it;
        L->nodes.splice(L->nodes.end(), std::move(nodes_cpy));

        int j = a(i, nodes.size()) + 1;
        int remaining = nodes.size();
        auto it = nodes.begin();
        Sublist* sublist = nullptr;
        std::list<T> sublist_elements;

        do {
            int prev_j = j;
            while (j >= 0 && 2 * A(i, j) > remaining) -- j;

            if (j != prev_j) {
                if (sublist != nullptr) {
                    sublist->elements = new List;
                    sublist->elements->sublist = sublist;
                    sublist->elements->id = id;
                    sublist->elements->i = i - 1;
                    sublist->elements->min_key = infinity;
                    sublist->elements->min_val = empty_val;

                    initialize_tail(sublist->elements, sublist_elements, id, i - 1);
                    sublist_elements.clear();
                    
                    L->tail.push_back(sublist);
                    sublist->sublist_it = std::prev(L->tail.end());
                }
                if (j == -1) {
                    sublist = nullptr;
                } else {
                    sublist = new Sublist;
                    sublist->list = L;
                    sublist->level = j;
                    sublist->head = false;
                }
            }    

            if (j == -1) {
                T singleton = *it;
                e[i][singleton] = -1;
                element_list[i][singleton] = L;
                L->tail_singletons.push_back(singleton);
                it ++;
            } else {
                int size = 2 * A(i, j);

                std::list<T> se_nodes;
                for (int j = 0; j < size; ++ j) {
                    se_nodes.push_back(*it);
                    it ++;
                }

                T superelement = se_nodes.back();
                sublist_elements.push_back(superelement);
                key[i-1][superelement] = infinity;
                val[i-1][superelement] = empty_val;

                for (auto n : se_nodes) {
                    e[i][n] = superelement;
                    update_key_at_level(superelement, i - 1, key[i][n], val[i][n]);
                }

                superelement_nodes[i][superelement] = std::move(se_nodes);

                remaining -= size;
            }
        } while (it != nodes.end());

        if (sublist != nullptr) {
            sublist->elements = new List;
            sublist->elements->sublist = sublist;
            sublist->elements->i = i - 1;
            sublist->elements->id = id;
            sublist->elements->min_key = infinity;
            sublist->elements->min_val = empty_val;

            initialize_tail(sublist->elements, sublist_elements, id, i - 1);
            
            L->tail.push_back(sublist);
            sublist->sublist_it = std::prev(L->tail.end());
        }

        for (auto n : nodes) update_key(L, key[i][n], val[i][n]);
    }

    List* find_list(T u, int i) {
        if (is_singleton(u, i)) {
            return element_list[i][u];
        } else {
            auto L = find_list(e[i][u], i - 1);
            return L->sublist->list;
        }
    }

    List* decrease_key(T u, K x, V v, int i) {
        update_key(key[i][u], val[i][u], x, v);
        if (is_singleton(u, i)) {
            auto L = element_list[i][u];
            update_key(L, x, v);
            return L;
        } else {
            auto eL = decrease_key(e[i][u], x, v, i - 1);
            auto L = eL->sublist->list;
            update_key(L, x, v);
            return L;
        }
    }

    std::pair<List*, List*> split(T x, H id1, H id2, int i) {
        if (is_singleton(x, i)) {
            List* L1 = element_list[i][x];
            L1->id = id1;
            List* L2 = new List;
            L2->id = id2;
            L2->i = i;
            L2->min_key = infinity;
            L2->min_val = empty_val;
            L2->nodes.splice(L2->nodes.begin(), L1->nodes, std::next(list_it[i][x]), L1->nodes.end());

            auto head_it = std::find(L1->head_singletons.begin(), L1->head_singletons.end(), x);
            auto tail_it = std::find(L1->tail_singletons.begin(), L1->tail_singletons.end(), x);

            if (head_it != L1->head_singletons.end()) {
                L2->head_singletons.splice(L2->head_singletons.end(), L1->head_singletons,
                                           std::next(head_it), L1->head_singletons.end());
                L2->head = std::move(L1->head);
                L2->tail = std::move(L1->tail);
                L2->tail_singletons = std::move(L1->tail_singletons);
            } else {
                assert(tail_it != L1->tail_singletons.end());

                L2->tail_singletons.splice(L2->tail_singletons.end(), L1->tail_singletons, 
                                           std::next(tail_it), L1->tail_singletons.end());
            }

            set_list_pointers(L2, i);

            calculate_min(L1, i);
            calculate_min(L2, i);

            return {L1, L2};
        } else {
            T superelement = e[i][x];
            auto [S1, S2] = split(superelement, id1, id2, i - 1);

            if (superelement == S1->nodes.front()) {
                S1->clear();
                calculate_min(S1, i - 1);
            } else {
                auto [S1p, S] = split(*std::prev(list_it[i-1][superelement]), id1, id2, i - 1);

                assert(S->nodes.size() == 1 && S->nodes.front() == superelement);
                assert(S1p == S1);

                delete S;
            }

            Sublist* S2sublist = new Sublist;
            S2sublist->elements = S2;
            S2->sublist = S2sublist;
            S2sublist->level = S1->sublist->level;
            S2sublist->head = S1->sublist->head;

            auto xit = std::find(superelement_nodes[i][superelement].begin(), 
                                 superelement_nodes[i][superelement].end(), x);
            std::list<T> nodes1(superelement_nodes[i][superelement].begin(), std::next(xit));
            std::list<T> nodes2(std::next(xit), superelement_nodes[i][superelement].end());

            auto L1 = S1->sublist->list;
            L1->id = id1;
            assert(L1->i == i);
            auto L2 = new List;
            L2->id = id2;
            L2->i = i;
            L2->min_key = infinity;
            L2->min_val = empty_val;

            L2->nodes.splice(L2->nodes.begin(), L1->nodes, 
                             std::next(list_it[i][superelement_nodes[i][superelement].back()]),
                             L1->nodes.end());
            
            L1->nodes.erase(list_it[i][superelement_nodes[i][superelement].front()], 
                            L1->nodes.end());

            if (S1->sublist->head) {
                auto split_it = S1->sublist->sublist_it;
                L2->head.splice(L2->head.end(), L1->head, std::next(split_it), L1->head.end());
                L2->head.push_front(S2sublist);
                S2sublist->sublist_it = L2->head.begin();
                S2sublist->head = true;
                L2->tail = std::move(L1->tail);
                L2->tail_singletons = std::move(L1->tail_singletons);
            } else {
                auto split_it = S1->sublist->sublist_it;
                L2->tail.splice(L2->tail.end(), L1->tail, std::next(split_it), L1->tail.end());
                L2->tail.push_front(S2sublist);
                S2sublist->sublist_it = L2->tail.begin();
                S2sublist->head = false;
                L2->tail_singletons = std::move(L1->tail_singletons);
            }

            set_list_pointers(L2, i);

            calculate_min(L1, i);
            calculate_min(L2, i);

            initialize_tail(L1, nodes1, id1, i);
            if (nodes2.size() > 0)
                initialize_head(L2, nodes2, id2, i);

            return {L1, L2};
        }
    }

    void calculate_min(List* L, int i) {
        L->min_key = infinity;
        L->min_val = empty_val;

        for (auto hs : L->head_singletons) update_key(L, key[i][hs], val[i][hs]);
        for (auto sl : L->head) update_key(L, sl->elements->min_key, sl->elements->min_val);
        for (auto sl : L->tail) update_key(L, sl->elements->min_key, sl->elements->min_val);
        for (auto ts : L->tail_singletons) update_key(L, key[i][ts], val[i][ts]);
    }

    void set_list_pointers(List* L, int i) {
        for (auto hs : L->head_singletons) element_list[i][hs] = L;
        for (auto sl : L->head) sl->list = L;
        for (auto sl : L->tail) sl->list = L;
        for (auto ts : L->tail_singletons) element_list[i][ts] = L;
    }

    bool is_singleton(T u, int i) {
        return e[i][u] == -1;
    }

    void update_key_at_level(T u, int i, K k, V v) {
        update_key(key[i][u], val[i][u], k, v);
    }

    void update_key(List* L, K key, V val) {
        update_key(L->min_key, L->min_val, key, val);
    }

    void update_key(K& old_key, V& old_val, K key, V val) {
        if (key < old_key) {
            old_key = key;
            old_val = val;
        }
    }

    void reset_keys(List* L) {
        for (auto n : L->nodes) {
            key[L->i][n] = infinity;
            val[L->i][n] = empty_val;
        }
        for (auto s : L->head) reset_keys(s->elements);
        for (auto s : L->tail) reset_keys(s->elements);
    }
};


template<typename T>
class ArrayPriorityQueue {
public:
    ArrayPriorityQueue() { reset(); }

    void scheduleEvent(int time, T event) {
        if (time < time) return;
        if (time >= event_queue.size())
            event_queue.resize(time * 2);
        event_queue[time].push_back(event);
    }

    T getEvent() {
        advanceToNext();
        return event_queue[time][it];
    }

    int timeNow() { return time; }

    void reset() {
        event_queue.clear();
        event_queue.resize(4);
        it = -1;
        time = 0;
    }

private:
    int time;
    size_t it;
    std::vector<std::vector<T>> event_queue;

    void advanceToNext() {
        it ++;
        while (it >= event_queue[time].size()) {
            time ++;
            it = 0;
        }
    }
};


class FenwickTree {
public:
    int sum(int index) {
        index ++;
        int res = 0;
        for (; index > 0; index -= index & (-index))
            res += T[index];
        return res;
    }

    void add(int index, int val) {
        index ++;
        for (; index < T.size(); index += index & (-index))
            T[index] += val;
    }

    void reset(int size) {
        T.resize(size + 1);
        for (int i = 0; i <= size; ++ i)
            T[i] = 0;
    }

private:
    std::vector<int> T;
};


} // namespace Koala