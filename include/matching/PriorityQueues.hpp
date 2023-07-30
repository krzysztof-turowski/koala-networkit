#pragma once

#include <set>
#include <vector>
#include <functional>
#include <iostream>
namespace Koala {

#define QUEUE_DEBUG 0

/**
 * Implementation of pq1 from the paper 'An O(EV log V) algorithm for finding a maximal 
 * weighted mathcing in general graphs' by Galil, Micali, Gabow. 
 * Used to store elements and their priorities and additionally allows for 
 * decreasing all priorities by a given amount efficiently.
 * 
 * @tparam E type of elmenents
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
        queue(cmp),
        modified_priority(size, 0) {}

    /**
     * Insert an element with given priority
     * 
     * @param value the value of inserted element
     * @param priority the priority of inserted element
    */
    void insert(E value, P priority) {
        // std::cerr << this << " insert element " << value << std::endl;
        modified_priority[value] = priority + Delta;
        queue.insert({ priority + Delta, value });
    }

    /**
     * Remove the element with provided value
     * 
     * @param value the value of element to remove
    */
    void remove(E value) { queue.erase(get_element(value)); }

    /**
     * Remove the minimum element
    */
    void remove_min() { queue.erase(queue.begin()); }

    /**
     * Check if queue is empty
     * 
     * @return 'true' if queue is empty, 'false' otherwise
    */
    bool empty() { return queue.empty(); }

    /**
     * Clear the queue
    */
    void clear() {
        queue.clear();
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
    */
    P current_priority(E value) {
        return modified_priority[value] - Delta;
    }

    /**
     * Find the current minimum element
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P>  find_min() {
        auto min = queue.begin();
        return { min->value, min->modified_priority - Delta };
    }

    /**
     * Call function for each element of queue. The function receives the value and priority of each element
     * 
     * @param handle the function to be called
    */
    void for_elements(const std::function<void(E,P)>& handle) {
        for (auto entry : queue) {
            handle(entry.value, entry.modified_priority - Delta);
        }
    }

    /**
     * Call function for each element of queue until it returns true or there are no more elements. The function receives the value and priority of each element
     * 
     * @param handle the function to be called
    */
    void for_elements_until(const std::function<bool(E,P)>& handle) {
        for (auto entry : queue) {
            if (handle(entry.value, entry.modified_priority - Delta)) {
                break;
            }
        }
    }

private:
    struct Entry { 
        P modified_priority; 
        E value; 
    };
    static bool cmp(const Entry& a, const Entry& b) { 
        return a.modified_priority < b.modified_priority 
            || (a.modified_priority == b.modified_priority && a.value < b.value); 
    }

    P Delta = 0;
    std::set<Entry, decltype(cmp)*> queue;
    std::vector<P> modified_priority;

    typename std::set<Entry>::iterator get_element(E value) {
        return queue.find({ modified_priority[value], value });
    }
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
    using ElementRef = Node*;

    H head;
    
    ConcatenableQueue(H head): head(head), root(nullptr) {}

    ConcatenableQueue(H head, E element, P priority): head(head), 
        root(new Node {
            this, nullptr, { nullptr, nullptr, nullptr },
            0, element, priority
        }) {}

    ConcatenableQueue(H head, Node* element): head(head) {
        // std::cerr << this << ": construct " << head << " "; element->print(); std::cerr << std::endl;
        if (element->height == 0) {
            root = new Node {
                this, nullptr, { element, nullptr, nullptr },
                1, element->min_element, element->min_priority
            };
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

    ConcatenableQueue(const ConcatenableQueue<H, E, P>& other) {
        // TODO make it copy contents or disable it
        root = other.root;
        head = other.head;
        fix_root();
        
        #if QUEUE_DEBUG
        check_consistency();
        #endif
    }

    void swap(ConcatenableQueue<H, E, P>& other) {
        std::swap(root, other.root);
        std::swap(head, other.head);
        fix_root();
        other.fix_root();
    }

    void fix_root() { if (root != nullptr) { root->parent = nullptr; root->queue = this; } }

    ConcatenableQueue& operator=(const ConcatenableQueue<H, E, P>& other) {
        if (this != &other) {
            auto temp(other);
            swap(temp);
        }
        return *this;
    }

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

    /**
     * Append element to the end of queue
     * 
     * @param element value of appended element
     * @param priority priority of appended element
     * @return reference to inserted element, used for spliting and finding queues. 
     *         Preserved by splits and concatenations.
    */
    ElementRef append(E element, P priority) {
        #if QUEUE_DEBUG
        std::cerr << this << " append " << element << std::endl;
        #endif

        if (empty()) {
            Node* new_node = new Node {
                this, nullptr, { nullptr, nullptr, nullptr },
                0, element, priority
            };    
            root = new Node {
                this, nullptr, { new_node, nullptr, nullptr },
                1, element, priority
            };
            new_node->parent = root;
            return new_node;
        }

        ElementRef last = last_element();
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
    ElementRef insert_after(ElementRef ref, E element, P priority) {
        #if QUEUE_DEBUG
        std::cerr << this << " insert after "; ref->print(); std::cerr << element << std::endl;
        #endif

        Node* new_node = new Node {
            this, nullptr, { nullptr, nullptr, nullptr },
            0, element, priority
        };

        add_child_after(ref->parent, ref, new_node);

        new_node->parent->update_min_until_root();

        #if QUEUE_DEBUG
        root->print(); std::cerr << std::endl;
        check_consistency();
        #endif

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
    ElementRef insert_before(ElementRef ref, E element, P priority) {
        #if QUEUE_DEBUG
        std::cerr << this << " insert " << element << " before "; ref->print(); 
        std::cerr << " in queue "; root->print(); std::cerr << std::endl;
        check_consistency();
        #endif

        Node* new_node = new Node {
            this, nullptr, { nullptr, nullptr, nullptr },
            0, element, priority
        };

        add_child_before(ref->parent, ref, new_node);

        new_node->parent->update_min_until_root();

        #if QUEUE_DEBUG
        root->print(); std::cerr << std::endl;
        check_consistency();
        #endif

        return new_node;
    }

    /**
     * Check if queue is empty
     * 
     * @return 'true' if queue is empty, 'false' otherwise
    */
    bool empty() { return root == nullptr; }

    /**
     * Remove an element
     * 
     * @param element_ref reference to the element to remove
    */
    void remove(ElementRef element_ref) {
        #if QUEUE_DEBUG
        std::cerr << this << " remove "; element_ref->print(); 
        std::cerr << " from "; root->print(); std::cerr << std::endl;
        #endif
        

        if (element_ref->parent == root && root->children_count() == 1) {
            delete root;
            delete element_ref;
            root = nullptr;
            return;
        }
        remove_node(element_ref);

        #if QUEUE_DEBUG
        root->print(); std::cerr << std::endl;
        check_consistency();
        #endif
    }

    /**
     * Find the current minimum element
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P> find_min() {
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

    void concat(ConcatenableQueue<H, E, P>&& other, H new_head, bool update_min = true) {
        #if QUEUE_DEBUG
        // std::cerr << "Concat: " << this << " and " << &other << "\n";
        // std::cerr << "  "; print_elements();
        // std::cerr << "  "; other.print_elements();
        #endif

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
            Node* new_root = new Node {
                this, nullptr, { left, right, nullptr },
                left->height + 1, left->min_element, left->min_priority
            };
            left->parent = new_root;
            right->parent = new_root;
            if (update_min) new_root->update_min();
            root = new_root;
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

        #if QUEUE_DEBUG
        // print_elements();
        check_consistency();
        #endif
    }

    std::pair<ConcatenableQueue<H, E, P>*, ConcatenableQueue<H, E, P>*> 
    split(ElementRef split_element, H head_left, H head_right) {
        #if QUEUE_DEBUG
        std::cerr << "Split " << this << " on "; 
        split_element->print(); 
        std::cerr << std::endl;
        print_elements();
        #endif

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

        #if QUEUE_DEBUG
        std::cerr << "left:  "; left->print_elements();
        std::cerr << "right: "; right->print_elements();
        left->check_consistency();
        right->check_consistency();
        #endif

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
            assert(false);
            // exit(1);
        }
        if (root->queue != this) {
            std::cerr << "root points at wrong queue\n";
            exit(1);
        }
        root->check_children();
    }

    struct Node {
        ConcatenableQueue* queue; // Correct in the root
        Node* parent;
        Node* children[3];
        int height;
        E min_element;
        P min_priority;

        int children_count() {
            if (children[0] == nullptr) return 0;
            if (children[1] == nullptr) return 1;
            if (children[2] == nullptr) return 2;
            return 3;
        }

        ConcatenableQueue* find_queue() {
            Node* iter = this;
            while (iter->parent != nullptr) iter = iter->parent;
            return iter->queue;
        }

        void update_min() {
            // std::cerr << "update min "; print();
            min_element = children[0]->min_element;
            min_priority = children[0]->min_priority;
            for (int i = 1; i < 3; ++ i) {
                if (children[i] == nullptr) break;
                if (children[i]->min_priority < min_priority) {
                    min_priority = children[i]->min_priority;
                    min_element = children[i]->min_element;
                }
            }
            // std::cerr << " -> " << min_element << std::endl;
        }

        void update_min_until_root() {
            // std::cerr << "updating min on path "; print(); std::cerr << std::endl;
            update_min();
            if (parent != nullptr) { 
                // std::cerr << "parent "; parent->print(); std::cerr << std::endl;
                parent->update_min_until_root();
            } 
            // else { std::cerr << "parent nullptr\n"; }
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
                out << children_count() /* << "|" << min_element */ << ":";
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
    };

private:
    Node* root;

    void insert_child(Node** children, int index, int child_count, Node* child) {
        for (int i = child_count; i > index; -- i)
            children[i] = children[i-1];
        children[index] = child;
    }

    void make_new_root(Node* a, Node* b) {
        // std::cerr << "make new root ("; 
        // a->print();
        // std::cerr << ", ";
        // b->print();
        // std::cerr << ")\n";

        Node* new_root = new Node { 
            this, nullptr, { a, b, nullptr },
            a->height + 1, a->min_element, a->min_priority
        };
        a->parent = new_root;
        b->parent = new_root;
        new_root->update_min();
        root = new_root;
    }

    void add_child_after(Node* target, Node* after, Node* child) {
        // std::cerr << "add child after ("; 
        // after->print();
        // std::cerr << ", ";
        // child->print();
        // std::cerr << ") in parent ";
        // if (target == nullptr) std::cerr << "nullptr";
        // else target->print();
        // std::cerr << "\n";

        if (target == nullptr) {
            make_new_root(after, child);
            return;
        }
        
        int after_index = target->child_index(after);
        insert_child_at_index(target, after_index + 1, child);
    }

    void add_child_before(Node* target, Node* before, Node* child) {
        // std::cerr << "add child before ("; 
        // before->print();
        // std::cerr << ", ";
        // child->print();
        // std::cerr << ") in parent ";
        // if (target == nullptr) std::cerr << "nullptr";
        // else target->print();
        // std::cerr << "\n";

        if (target == nullptr) {
            make_new_root(child, before);
            return;
        }
        
        int before_index = target->child_index(before);
        insert_child_at_index(target, before_index, child);
    }

    void insert_child_at_index(Node* target, int index, Node* child) {
        // std::cerr << "add child "; 
        // child->print();
        // std::cerr << " at index " << index <<  " in target ";
        // if (target == nullptr) std::cerr << "nullptr";
        // else target->print();
        // std::cerr << "\n";

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
            Node* new_node = new Node {
                this, target->parent, { four_children[2], four_children[3], nullptr},
                target->height, four_children[2]->min_element, four_children[2]->min_priority
            };
            four_children[2]->parent = new_node;
            four_children[3]->parent = new_node;
            target->update_min();
            new_node->update_min();
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

    ElementRef first_element() {
        Node* iter = root;
        while (!iter->is_leaf()) {
            iter = iter->leftmost_child();
        }
        return iter;
    }

    ElementRef last_element() {
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
 * @tparam E type of elmenents
 * @tparam P type of priority values
*/
template<typename E, typename P>
class PriorityQueue2 {
public:
    struct Group { 
        Group(bool active, P Delta_last, P Delta_group): 
            elements(this), active(active), Delta_last(Delta_last), Delta_group(Delta_group) {}
        
        Group(ConcatenableQueue<Group*, E, P>&& elements_, 
                bool active, P Delta_last, P Delta_group): 
            elements(std::move(elements_)),
            active(active), Delta_last(Delta_last), Delta_group(Delta_group) {

            elements.head = this;
        }

        ConcatenableQueue<Group*, E, P> elements;
        bool active;
        P Delta_last;
        P Delta_group;

        bool empty() {
            return elements.empty() || elements.find_min().second == dummy_element_priority;
        }

        std::pair<E, P> find_min() {
            auto [value, priority] = elements.find_min();
            return {value, priority - Delta_group};
        }

        void update_Delta(P Delta) {
            if (active) {
                Delta_group = Delta_group + Delta - Delta_last;
            }
            Delta_last = Delta;
        }
    };
    
    /**
     * Creates a queue that stores elements between 0 and size - 1.
     * 
     * @param size the maximum value of elements that can be stored
    */
    PriorityQueue2(E size): group_minima(size), elements(size, nullptr) {}

    /**
     * Append a dummy element to the end group
     * 
     * @param value the value of inserted dummy element
     * @param group the group to which the new element belongs
    */
    void append_dummy(E value, Group* group) {
        elements[value] = group->elements.append(value, dummy_element_priority);

        #if QUEUE_DEBUG
        std::cerr << "=========================================\n";
        std::cerr << "append dummy " << value << " to group ";
        group->elements.print_elements();
        group->elements.check_consistency();
        #endif
    }

    /**
     * Append an element with given priority to the end group
     * 
     * @param value the value of inserted element
     * @param priority the priority of inserted element
     * @param group the group to which the new element belongs
    */
    void append(E value, P priority, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group->find_min();

        #if QUEUE_DEBUG
        std::cerr << "=========================================\n";
        std::cerr << "append " << value << " to group " << group->active << " ";
        group->elements.print_elements();
        std::cerr << "with min " << last_min << " " << last_min_priority << std::endl;
        #endif

        elements[value] = group->elements.append(value, priority + group->Delta_group);
        
        #if QUEUE_DEBUG
        group->elements.check_consistency();
        #endif

        // Check if the new element is the minimum in group
        if (group->active && group->find_min().first == value) {
            if (last_min_priority != dummy_element_priority)
                group_minima.remove(last_min);
            group_minima.insert(value, priority);
        }
    }

    /**
     * Insert an element with given priority into group before element
     * 
     * @param value the value of inserted element
     * @param priority the priority of inserted element
     * @param before the value of element before which the new element is inserted
     * @param group the group to which the new element belongs
    */
    void insert_before(E value, P priority, E before, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group->find_min();

        #if QUEUE_DEBUG
        std::cerr << "=========================================\n";
        std::cerr << "insert " << value << " with priority " << priority
                  << " before " << before << " in group " << group->active << " ";
        group->elements.print_elements();
        std::cerr << "with min " << last_min << " " << last_min_priority << std::endl;
        group_minima.for_elements([] (E u, P w) {
            std::cerr << u << " " << w << std::endl;
        });
        #endif

        elements[value] = group->elements.insert_before(elements[before], value, priority + group->Delta_group);

        #if QUEUE_DEBUG
        group->elements.check_consistency();
        #endif
        
        // Check if the new element is the minimum in group
        if (group->active && group->find_min().first == value) {
            if (last_min_priority != dummy_element_priority)
                group_minima.remove(last_min);
            group_minima.insert(value, priority);
        }
    }

    /**
     * Remove an element
     * 
     * @param value the value of the element to remove
    */
    void remove(E value) { 
        Group* group = elements[value]->find_queue()->head;
        group->update_Delta(Delta);
        auto [last_min, last_min_priority] = group->find_min();

        #if QUEUE_DEBUG
        std::cerr << "=========================================\n";
        std::cerr << "remove " << value << " from group ";
        group->elements.print_elements();
        std::cerr << "with min " << last_min << " " << last_min_priority << std::endl;
        #endif
        
        group->elements.remove(elements[value]);

        #if QUEUE_DEBUG
        group->elements.check_consistency();
        #endif

        // Check if the removed element was the minimum in group
        if (last_min == value && group->active) {
            group_minima.remove(value);
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
    bool has_active_elements() {
        return !group_minima.empty();
    }

    /**
     * Find the current minimum element in any active group
     * 
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<E, P>  find_min() {
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
        #if QUEUE_DEBUG
        auto [min, min_p] = group->find_min();
        std::cerr << "=========================================\n";
        std::cerr << "delete group " << group->active << " ";
        group->elements.print_elements();
        std::cerr << "with min " << min << " " << min_p << std::endl;
        group_minima.for_elements([] (E u, P w) {
            std::cerr << u << " " << w << std::endl;
        });
        #endif

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

        #if QUEUE_DEBUG
        auto [min, min_p] = group->find_min();
        std::cerr << "=========================================\n";
        std::cerr << "change status to " << active << " of group ";
        group->elements.print_elements();
        std::cerr << "with min " << min << " " << min_p << std::endl;
        group_minima.for_elements([] (E u, P w) {
            std::cerr << u << " " << w << std::endl;
        });
        #endif

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
     * @param value the value of the element according to which the groups are to be split
     * @return a pair containing the two new groups, the first contains elements from beginning of group
     *      until and containing the split element, the second one contains elements after the element
    */
    std::pair<Group*, Group*> split_group(Group* group, E value) {
        group->update_Delta(Delta);

        #if QUEUE_DEBUG
        auto [min, min_p] = group->find_min();
        std::cerr << "=========================================\n";
        std::cerr << "split at " << value << " the gruop " << group->active << " ";
        group->elements.print_elements();
        std::cerr << "with min " << min << " " << min_p << std::endl;
        group_minima.for_elements([] (E u, P w) {
            std::cerr << u << " " << w << std::endl;
        });
        #endif

        if (group->active && !group->empty()) {
            auto [v, p] = group->find_min();
            group_minima.remove(v);
        }

        auto [elements_left, elements_right] = group->elements.split(elements[value], group, group);

        #if QUEUE_DEBUG
        elements_left->print_elements();
        elements_right->print_elements();
        #endif

        Group* group_left  = new Group(std::move(*elements_left),
                                group->active, group->Delta_last, group->Delta_group);
        Group* group_right = new Group(std::move(*elements_right),
                                group->active, group->Delta_last, group->Delta_group);
        
        #if QUEUE_DEBUG
        std::cerr << group_left << " " << group_right << std::endl;
        group_left->elements.print_elements();
        group_right->elements.print_elements();
        #endif

        if (group->active) {
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

    void shift_group(Group* group, E value) {
        group->update_Delta(Delta);

        auto [elements_left, elements_right] = group->elements.split(elements[value], group, group);

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
    static constexpr P dummy_element_priority = std::numeric_limits<P>::max();

    PriorityQueue1<E, P> group_minima;
    std::vector<typename ConcatenableQueue<Group*, E, P>::ElementRef> elements;
    
    P Delta = 0;
};


// template<typename E, typename P>
// class FHeap {
// public:
//     void insert(E val, P priority) {
//         nodes.push_back({ val, priority });
//     }

//     P find_min() {

//     }

//     void delete_min() {
//         if (min_node == nodes.end()) return;
//         nodes.splice(nodes.end(), min_node->children);
//         nodes.erase(min_node);
//         consolidaate();
//     }

//     void decrease_cost(E val, P priority) {

//     }

//     void merge(FHeap& other) {
//         nodes.splice(nodes.end(), other.nodes);
//         min_node = min_node->priority < other.min_node->priority ? 
//             min_node : other.min_node;
//     }

// private:
//     struct Node {
//         Node(): value(0), priority(0) {}
//         Node(E val, P priority): value(val), priority(priority) {}

//         E value;
//         P priority;
//         std::list<Node> children;

//         int degree() { return children.size(); }
//     };

//     std::list<Node> nodes;
//     std::list<Node>::iterator min_node;

//     int max_degree() {
//         int res = 0;
//         for (auto n : nodes)
//     }

//     void consolidate() {
//         int D = max_degree();
//         Node nodes[D + 1];

//         for (auto node : nodes) {
//             Node new_node = node;
//             int deg = new_node.degree();

//             while (nodes[])
//         }

//         nodes.clear();
//     }
// };


// template<typename E, typename P>
// class SplittingList {
// public:
//     void insert(E val, P priority) {
//         nodes.push_back({ val, priority });
//     }

//     P find_min() {

//     }

//     void delete_min() {
//         if (min_node == nodes.end()) return;
//         nodes.splice(nodes.end(), min_node->children);
//         nodes.erase(min_node);
//         consolidaate();
//     }

//     void decrease_cost(E val, P priority) {

//     }

//     void merge(FHeap& other) {
//         nodes.splice(nodes.end(), other.nodes);
//         min_node = min_node->priority < other.min_node->priority ? 
//             min_node : other.min_node;
//     }

// private:
//     struct Node {
//         Node(): value(0), priority(0) {}
//         Node(E val, P priority): value(val), priority(priority) {}

//         E value;
//         P priority;
//         std::list<Node> children;

//         int degree() { return children.size(); }
//     };

//     std::list<Node> nodes;
//     std::list<Node>::iterator min_node;

//     int max_degree() {
//         int res = 0;
//         for (auto n : nodes)
//     }

//     void consolidate() {
//         int D = max_degree();
//         Node nodes[D + 1];

//         for (auto node : nodes) {
//             Node new_node = node;
//             int deg = new_node.degree();

//             while (nodes[])
//         }

//         nodes.clear();
//     }
// };

template<typename N>
class UnionFind {
public:
    void reset() {
        name.clear();
        parent.clear();
    }

    int create(N x_name) {
        int x = name.size();
        name.push_back(x_name);
        parent.push_back(x);
        return x;
    }

    int link(int x, N y_name) {
        int y = name.size();
        name.push_back(y_name);
        parent.push_back(y);
        return y;
    }

    N find(int x) {
        return name[find_root(x)];
    }

private:
    std::vector<N> name;
    std::vector<int> parent;

    int find_root(int x) {
        return parent[x] == x ? x : (parent[x] = find_root(parent[x]));
    }
};


} // namespace Koala