#pragma once

namespace Koala {

/**
 * Implementation of concatenable queues. Stores elements with their priorities.
 * Preserves the order of elements that is independent from their priorities.
 * Allows appending and inserting elements as well as concatenating queues and spliting them.
 * Each queue has an associated id;
 *
 * @tparam Id type of the id associated with each queue
 * @tparam Element type of elements
 * @tparam Priority type of priority values
*/
template<class Element, class Priority, class Id>
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
    handle_type append(Element element, Priority priority) {
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
    handle_type insert_after(handle_type ref, Element element, Priority priority) {
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
    handle_type insert_before(handle_type ref, Element element, Priority priority) {
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
    void decrease_priority(handle_type ref, Priority priority) {
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
    std::pair<Element, Priority> find_min() const {
        return { root->min_element, root->min_priority };
    }

    /**
     * Concatenate two queues
     *
     * @param left the left queue to concatenate
     * @param right the right queue to concatenat
     * @return pointer to a new queue that is a concatenation of the two provided queues
    */
    static ConcatenableQueue* concat(
            ConcatenableQueue&& left,
            ConcatenableQueue&& right,
            Id head, bool update_min = true) {
        auto new_left = new ConcatenableQueue(std::move(left));
        new_left->concat(std::move(right), head, update_min);
        return new_left;
    }

    /**
     * Concatenate the provieded queue to the end of the current one while updating the head value.
     *
     * @param other the left queue to concatenate
     * @param new_head the head value of the new queue
    */
    void concat(ConcatenableQueue&& other, Id new_head, bool update_min = true) {
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

    void for_each(const std::function<void(Element, Priority)>& handle) {
        if (root != nullptr) root->for_each(handle);
    }

    void set_all_priorities(Priority priority) {
        root->set_all_priorities(priority);
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
        Element min_element;
        Priority min_priority;

        Node(ConcatenableQueue* queue, Element element, Priority priority):
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

        void for_each(const std::function<void(Element, Priority)>& handle) {
            if (is_leaf()) {
                handle(min_element, min_priority);
                return;
            }
            for (int i = 0; i < 3; ++i) {
                if (children[i] == nullptr) return;
                children[i]->for_each(handle);
            }
        }

        void set_all_priorities(Priority priority) {
            min_priority = priority;
            for (int i = 0; i < 3; ++i) {
                if (children[i] == nullptr) break;
                children[i]->set_all_priorities(priority);
            }
            if (children[0] != nullptr) {
                min_element = children[0]->min_element;
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

    ConcatenableQueue(Id id, Element element, Priority priority):
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

} /* namespace Koala */
