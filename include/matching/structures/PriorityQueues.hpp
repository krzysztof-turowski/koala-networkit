#pragma once

#include <set>
#include <list>
#include <vector>
#include <functional>

#include <matching/structures/HeapWithRemove.hpp>

namespace Koala {

/**
 * Implementation of pq1 from the paper 'An O(EV log V) algorithm for finding a maximal
 * weighted mathcing in general graphs' by Galil, Micali, Gabow.
 * Used to store elements and their priorities and additionally allows for
 * decreasing all priorities by a given amount efficiently.
 *
 * @tparam E type of elements
 * @tparam V type of values
 * @tparam P type of priority values
*/
template<typename E, typename V, typename P>
class PriorityQueue1 {
 public:
    /**
     * Creates a queue that stores elements between 0 and size - 1.
     *
     * @param size the maximum value of elements that can be stored
    */
    explicit PriorityQueue1(E size):
        heap(),
        element_handle(size),
        in_heap(size, false),
        modified_priority(size, 0),
        element_value(size) {}

    /**
     * Insert an element with given priority
     *
     * @param element the value of inserted element
     * @param priority the priority of inserted element
    */
    void insert(E element, V value, P priority) {
        if (!in_heap[element]) {
            modified_priority[element] = priority + Delta;
            element_handle[element] = heap.push({ priority + Delta, element });
            in_heap[element] = true;
            element_value[element] = value;
        } else if (priority < current_priority(element)) {
            remove(element);
            modified_priority[element] = priority + Delta;
            element_handle[element] = heap.push({ priority + Delta, element });
            element_value[element] = value;
        }
    }

    /**
     * Remove the element with provided value
     *
     * @param element the value of element to remove
    */
    void remove(E element) {
        heap.erase(element_handle[element]);
        in_heap[element] = false;
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
     * @param element the value of element to check
     * @returns current priority of the provided element
    */
    P current_priority(E element) const {
        return modified_priority[element] - Delta;
    }

    /**
     * Find the current minimum element
     *
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<V, P> find_min() const {
        auto [priority, element] = heap.top();
        return { element_value[element], priority - Delta };
    }

    /**
     * Find the current minimum element
     *
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::tuple<E, V, P> find_min_element() const {
        auto [priority, element] = heap.top();
        return { element, element_value[element], priority - Delta };
    }

    /**
     * Call function for each element of queue. The function receives each element with it's
     * priority and the associated value
     *
     * @param handle the function to be called
    */
    void for_elements(const std::function<void(E, V, P)>& handle) {
        for (auto entry : heap) {
            handle(entry.second, element_value[entry.second], entry.first - Delta);
        }
    }

 private:
    P Delta = 0;
    using pair_heap = HeapWithRemove<std::pair<P, E>>;
    pair_heap heap;
    std::vector<typename pair_heap::handle_type> element_handle;
    std::vector<bool> in_heap;
    std::vector<P> modified_priority;
    std::vector<V> element_value;
};

/**
 * Implementation of concatenable queues. Stores elements with their priorities.
 * Preserves the order of elements that is independent from their priorities.
 * Allows appending and inserting elements as well as concatenating queues and spliting them.
 * Each queue has an associated id;
 *
 * @tparam Id type of the id value stored in each queue
 * @tparam E type of element values
 * @tparam P type of priority values
*/
template<class Id, class E, class P>
class ConcatenableQueue {
 public:
    struct Node;
    using handle_type = Node*;

    Id root_id;

    explicit ConcatenableQueue(Id head): root_id(head), root(nullptr) {}

    ConcatenableQueue(ConcatenableQueue&& other) {
        root = other.root;
        root_id = other.root_id;
        other.root = nullptr;
        fix_root();
    }

    ConcatenableQueue& operator=(ConcatenableQueue&& other) {
        if (this != &other) {
            swap(other);
        }
        return *this;
    }

    ~ConcatenableQueue() {
        if (root != nullptr) {
            root->delete_children();
            delete root;
        }
    }

    // Delete copy constructor as we don't want it used
    // TODO(rkilar): implement copying maybe
    ConcatenableQueue(const ConcatenableQueue& other) = delete;
    ConcatenableQueue& operator=(const ConcatenableQueue& other) = delete;

    void swap(ConcatenableQueue& other) {
        std::swap(root, other.root);
        std::swap(root_id, other.root_id);
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
     * Decreases priority of an element to a new value if it's better than the current one.
     *
     * @param ref reference to the element whose priority is to be updated
     * @param priority new priority
    */
    void decrease_priority(handle_type ref, P priority) {
        if (priority < ref->min_priority) {
            ref->min_priority = priority;
            if (ref->parent != nullptr) {
                ref->parent->update_min_until_root();
            }
        }
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

    /**
     * Concatenate two queues
     *
     * @param left the left queue to concatenate
     * @param right the right queue to concatenat
     * @return pointer to a new queue that is a concatenation of the two provided queues
    */
    static ConcatenableQueue<Id, E, P>* concat(
            ConcatenableQueue<Id, E, P>&& left, ConcatenableQueue<Id, E, P>&& right,
            Id head, bool update_min = true) {
        auto new_left = new ConcatenableQueue<Id, E, P>(std::move(left));
        new_left->concat(std::move(right), head, update_min);
        return new_left;
    }

    /**
     * Concatenate the provieded queue to the end of the current one while updating the head value.
     *
     * @param other the left queue to concatenate
     * @param new_head the head value of the new queue
    */
    void concat(ConcatenableQueue<Id, E, P>&& other, Id new_head, bool update_min = true) {
        if (other.empty()) { root_id = new_head; return; }
        if (empty()) {
            *this = std::move(other);
            root_id = new_head;
            return;
        }
        root_id = new_head;
        Node* left = root;
        cut(left);
        Node* right = other.root;
        cut(right);
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

    std::pair<ConcatenableQueue*, ConcatenableQueue*>
    split(handle_type split_element, Id id_left, Id id_right) {
        Node* prev = split_element;
        Node* iter = split_element->parent;
        ConcatenableQueue* left = new ConcatenableQueue(id_left, split_element);
        ConcatenableQueue* right = new ConcatenableQueue(id_right);

        while (iter != nullptr) {
            int child_count = iter->children_count();
            int path = iter->child_index(prev);

            for (int i = path - 1; i >= 0; --i) {
                auto l = concat(iter->children[i], std::move(*left), id_left, false);
                delete left;
                left = l;
            }

            for (int i = path + 1; i < child_count; ++i) {
                auto r = concat(std::move(*right), iter->children[i], id_right, false);
                delete right;
                right = r;
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

    struct Node {
        ConcatenableQueue* find_queue() {
            Node* iter = this;
            while (iter->parent != nullptr) iter = iter->parent;
            return iter->queue;
        }

     private:
        // Pointer to the queue which is correct in the root
        ConcatenableQueue* queue;
        Node* parent;
        Node* children[3];
        int height;
        E min_element;
        P min_priority;

        Node(ConcatenableQueue* queue, E element, P priority):
            queue(queue),
            parent(nullptr),
            children{nullptr, nullptr, nullptr},
            height(0),
            min_element(element),
            min_priority(priority) {}

        static Node* from_child(ConcatenableQueue* queue, Node* child) {
            Node* node = new Node(queue, child->min_element, child->min_priority);
            node->children[0] = child;
            node->height = child->height + 1;
            child->parent = node;
            return node;
        }

        static Node* from_children(ConcatenableQueue* queue, Node* left, Node* right) {
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
            for (int i = 1; i < 3; ++i) {
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
            for (int i = 0; i < 3; ++i)
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

        void for_each(const std::function<void(E, P)>& handle) {
            if (is_leaf()) {
                handle(min_element, min_priority);
                return;
            }
            for (int i = 0; i < 3; ++i) {
                if (children[i] == nullptr) return;
                children[i]->for_each(handle);
            }
        }

        void delete_children() {
            for (int i = 0; i < 3; ++i)
                if (children[i] != nullptr) {
                    children[i]->delete_children();
                    delete children[i];
                }
        }

        friend class ConcatenableQueue;
    };

 private:
    Node* root;

    ConcatenableQueue(Id id, E element, P priority):
        root_id(id),
        root(new Node(element, priority)) {}

    ConcatenableQueue(Id id, Node* element):
            root_id(id) {
        if (element->height == 0) {
            root = Node::from_child(this, element);
            element->parent = root;
        } else {
            root = element;
            root->parent = nullptr;
        }

        fix_root();
    }

    static ConcatenableQueue* concat(
            ConcatenableQueue&& left, ConcatenableQueue::Node* right,
            Id id, bool update_min = true) {
        ConcatenableQueue right_que(id, right);
        return concat(std::move(left), std::move(right_que), id, update_min);
    }

    static ConcatenableQueue* concat(
                ConcatenableQueue::Node* left, ConcatenableQueue&& right,
            Id id, bool update_min = true) {
        ConcatenableQueue left_que(id, left);
        return concat(std::move(left_que), std::move(right), id, update_min);
    }

    void fix_root() {
        if (root != nullptr) {
            root->parent = nullptr;
            root->queue = this;
        }
    }

    void insert_child(Node** children, int index, int child_count, Node* child) {
        for (int i = child_count; i > index; --i)
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
        for (int i = child_index; i < children - 1; ++i)
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
        bool has_elements() const {
            return !elements.empty();
        }

        std::tuple<E, V, P> find_min() const {
            auto [element, pv] = elements.find_min();
            auto [priority, value] = pv;
            return {element, value, priority - Delta_group};
        }

        bool is_active() const {
            return active;
        }

     private:
        using ElementQueue = ConcatenableQueue<Group*, E, std::pair<P, V>>;

        Group(bool active, P Delta_last, P Delta_group):
            elements(this), active(active), Delta_last(Delta_last), Delta_group(Delta_group) {}

        Group(ElementQueue&& elements_, bool active, P Delta_last, P Delta_group):
                elements(std::move(elements_)),
                active(active),
                Delta_last(Delta_last),
                Delta_group(Delta_group) {
            elements.root_id = this;
        }

        void update_Delta(P Delta) {
            if (active) {
                Delta_group = Delta_group + Delta - Delta_last;
            }
            Delta_last = Delta;
        }

        ConcatenableQueue<Group*, E, typename std::pair<P, V>> elements;
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
    void append(E element, V value, P priority, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_value, last_min_priority] = group_minimum(group);

        elements[element] = group->elements.append(element, {priority + group->Delta_group, value});

        // Check if the new element is the minimum in group
        if (group->active && std::get<0>(group->find_min()) == element && value != no_value) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(element, value, priority);
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
    void insert_before(E element, V value, P priority, E before, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_value, last_min_priority] = group_minimum(group);

        elements[element] = group->elements.insert_before(elements[before], element,
            {priority + group->Delta_group, value});

        // Check if the new element is the minimum in group
        if (group->active && std::get<0>(group->find_min()) == element && value != no_value) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(element, value, priority);
        }
    }

    /**
     * Change priority of an element to a new value if it's lower than the current one, in which
     * case also update the associated value
     *
     * @param element the value of the element whose priority is to be changed
     * @param value value associated with the new priority
     * @param priority the new priority of the element
     * @param group the group to which the element belongs
    */
    void decrease_priority(E element, V value, P priority, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_value, last_min_priority] = group_minimum(group);

        if (last_min == element && last_min_priority <= priority) return;

        group->elements.decrease_priority(
            elements[element], {priority + group->Delta_group, value});

        // Check if the new element is the minimum in group
        if (group->active && std::get<0>(group->find_min()) == element && value != no_value) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(element, value, priority);
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
     * @return a pair containing the value associated with the lowest priority and the priority
    */
    std::pair<V, P> find_min() const {
        return group_minima.find_min();
    }

    /**
     * Find the current minimum element in any active group
     *
     * @return a tuple containing the element with the lowest priority, the associated value and
     * the priority
    */
    std::tuple<E, V, P> find_min_element() const {
        return group_minima.find_min_element();
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
        if (group->active && !is_empty(group)) {
            auto [e, v, p] = group->find_min();
            group_minima.remove(e);
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

        if (active && !is_empty(group)) {
            auto [e, v, p] = group->find_min();
            group_minima.insert(e, v, p);
        } else if (!active && !is_empty(group)) {
            auto [e, v, p] = group->find_min();
            group_minima.remove(e);
        }
    }

    /**
     * Split a given group acording to the given element
     *
     * @param group the group to be split
     * @param element the value of the element according to which the groups are to be split
     * @return a pair containing the two new groups, the first contains elements from beginning
     *      of group until and containing the split element, the second one contains elements after
     *      the element
    */
    std::pair<Group*, Group*> split_group(Group* group, E element) {
        group->update_Delta(Delta);

        if (group->active && !is_empty(group)) {
            auto [e, v, p] = group->find_min();
            group_minima.remove(e);
        }

        auto [left, right] = group->elements.split(elements[element], group, group);

        Group* group_left  = new Group(std::move(*left),
                                group->active, group->Delta_last, group->Delta_group);
        Group* group_right = new Group(std::move(*right),
                                group->active, group->Delta_last, group->Delta_group);

        delete left;
        delete right;

        if (group->active) {
            if (!is_empty(group)) {
                auto [e, v, p] = group->find_min();
                group_minima.remove(e);
            }
            if (!is_empty(group_left)) {
                auto [e, v, p] = group_left->find_min();
                group_minima.insert(e, v, p);
            }
            if (!is_empty(group_right)) {
                auto [e, v, p] = group_right->find_min();
                group_minima.insert(e, v, p);
            }
        }
        delete group;
        return {group_left, group_right};
    }


    bool is_empty(Group* group) {
        return !group->has_elements() || std::get<1>(group->find_min()) == no_value;
    }

    void shift_group(Group* group, E element) {
        group->update_Delta(Delta);

        auto [left, right] = group->elements.split(elements[element], group, group);

        right->concat(std::move(*left), group);
        group->elements = std::move(*right);
        delete left;
        delete right;
    }

    void for_each_in_group(Group* group, const std::function<void(E, V, P)>& handle) {
        group->update_Delta(Delta);
        group->elements.for_each([&handle, group] (E e, std::pair<P, V> pv) {
            handle(e, pv.second, pv.first - group->Delta_group);
        });
    }

    void for_group_minima(const std::function<void(E, V, P)>& handle) {
        group_minima.for_elements(handle);
    }

 private:
    std::tuple<E, V, P> group_minimum(Group* group) {
        return is_empty(group) ?
            std::make_tuple(no_element, no_value, infinite_priority) : group->find_min();
    }

    P Delta = 0;
    E no_element;
    V no_value;
    P infinite_priority;
    PriorityQueue1<E, V, P> group_minima;
    std::vector<typename Group::ElementQueue::handle_type> elements;
};

} /* namespace Koala */
