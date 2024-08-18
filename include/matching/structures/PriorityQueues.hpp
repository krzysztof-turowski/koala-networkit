#pragma once

#include <set>
#include <list>
#include <vector>
#include <functional>

#include <matching/structures/ConcatenableQueue.hpp>
#include <matching/structures/HeapWithRemove.hpp>

namespace Koala {

/**
 * Implementation of pq1 from the paper 'An O(nm log n) algorithm for finding a maximal
 * weighted mathcing in general graphs' by Galil, Micali, Gabow.
 * Used to store elements and their priorities and additionally allows for
 * decreasing all priorities by a given amount efficiently.
 *
 * @tparam Element type of elements
 * @tparam Priority type of priority values
 * @tparam Value type of values
*/
template<typename Element, typename Priority, typename Value>
class PriorityQueue1 {
 public:
    /**
     * Creates a queue that stores elements between 0 and size - 1.
     *
     * @param size the maximum value of elements that can be stored
    */
    explicit PriorityQueue1(Element size):
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
    void insert(Element element, Value value, Priority priority) {
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
    void remove(Element element) {
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
    void decrease_all_priorities(Priority delta) { Delta += delta; }

    /**
     * Check the current priority of element
     *
     * @param element the value of element to check
     * @returns current priority of the provided element
    */
    Priority current_priority(Element element) const {
        return modified_priority[element] - Delta;
    }

    /**
     * Find the current minimum element
     *
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::pair<Value, Priority> find_min() const {
        auto [priority, element] = heap.top();
        return { element_value[element], priority - Delta };
    }

    /**
     * Find the current minimum element
     *
     * @return a pair containing the value of an element with lowest priority it's priority
    */
    std::tuple<Element, Value, Priority> find_min_element() const {
        auto [priority, element] = heap.top();
        return { element, element_value[element], priority - Delta };
    }

    /**
     * Call function for each element of queue. The function receives each element with it's
     * priority and the associated value
     *
     * @param handle the function to be called
    */
    void for_elements(const std::function<void(Element, Value, Priority)>& handle) {
        for (auto entry : heap) {
            handle(entry.second, element_value[entry.second], entry.first - Delta);
        }
    }

 private:
    Priority Delta = 0;
    using pair_heap = HeapWithRemove<std::pair<Priority, Element>>;
    pair_heap heap;
    std::vector<typename pair_heap::handle_type> element_handle;
    std::vector<bool> in_heap;
    std::vector<Priority> modified_priority;
    std::vector<Value> element_value;
};

/**
 * Implementation of pq2 from the paper 'An O(nm log n) algorithm for finding a maximal
 * weighted mathcing in general graphs' by Galil, Micali, Gabow.
 * Used to store elements and their priorities.
 * The elements are divided into groups that can be either active or nonactive.
 * Allows for changing the status of a group, changing priorities off all elements in active groups,
 * splitting groups
 *
 * @tparam Element type of elements
 * @tparam Value type of values
 * @tparam Priority type of priority values
*/
template<typename Element, typename Priority, typename Value>
class PriorityQueue2 {
 public:
    class Group {
     public:
        bool has_elements() const {
            return !elements.empty();
        }

        std::tuple<Element, Value, Priority> find_min() const {
            auto [element, pv] = elements.find_min();
            auto [priority, value] = pv;
            return {element, value, priority - Delta_group};
        }

        bool is_active() const {
            return active;
        }

     private:
        using ElementlementQueue = ConcatenableQueue<
            Element, typename std::pair<Priority, Value>, Group*>;

        Group(bool active, Priority Delta_last, Priority Delta_group):
            elements(this), active(active), Delta_last(Delta_last), Delta_group(Delta_group) {}

        Group(ElementlementQueue&& elements_, bool active, Priority Delta_last,
            Priority Delta_group):
                elements(std::move(elements_)),
                active(active),
                Delta_last(Delta_last),
                Delta_group(Delta_group) {
            elements.root_id = this;
        }

        void update_Delta(Priority Delta) {
            if (active) {
                Delta_group = Delta_group + Delta - Delta_last;
            }
            Delta_last = Delta;
        }

        ElementlementQueue elements;
        bool active;
        Priority Delta_last;
        Priority Delta_group;

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
    PriorityQueue2(Element size, Element no_element, Value no_value, Priority infinite_priority):
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
    void append(Element element, Value value, Priority priority, Group* group) {
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
    void insert_before(Element e, Value value, Priority priority, Element before, Group* group) {
        group->update_Delta(Delta);
        auto [last_min, last_value, last_min_priority] = group_minimum(group);

        elements[ecvt_r] = group->elements.insert_before(elements[before], e,
            {priority + group->Delta_group, value});

        // Check if the new element is the minimum in group
        if (group->active && std::get<0>(group->find_min()) == e && value != no_value) {
            if (last_min != no_element)
                group_minima.remove(last_min);
            group_minima.insert(e, value, priority);
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
    void decrease_priority(Element element, Value value, Priority priority, Group* group) {
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
    std::pair<Value, Priority> find_min() const {
        return group_minima.find_min();
    }

    /**
     * Find the current minimum element in any active group
     *
     * @return a tuple containing the element with the lowest priority, the associated value and
     * the priority
    */
    std::tuple<Element, Value, Priority> find_min_element() const {
        return group_minima.find_min_element();
    }

    /**
     * Decrease all priorities of current elements in active groups by a provided value
     *
     * @param delta the amount by which to decrease priorities
    */
    void decrease_all_priorities(Priority delta) {
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
    std::pair<Group*, Group*> split_group(Group* group, Element element) {
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

    void shift_group(Group* group, Element element) {
        group->update_Delta(Delta);

        auto [left, right] = group->elements.split(elements[element], group, group);

        right->concat(std::move(*left), group);
        group->elements = std::move(*right);
        delete left;
        delete right;
    }

    using ElementHandlerFunction = std::function<void(Element, Value, Priority)>;

    void for_each_in_group(Group* group, const ElementHandlerFunction& handle) {
        group->update_Delta(Delta);
        group->elements.for_each([&handle, group] (Element e, std::pair<Priority, Value> pv) {
            handle(e, pv.second, pv.first - group->Delta_group);
        });
    }

    void for_group_minima(const ElementHandlerFunction& handle) {
        group_minima.for_elements(handle);
    }

 private:
    std::tuple<Element, Value, Priority> group_minimum(Group* group) {
        return is_empty(group) ?
            std::make_tuple(no_element, no_value, infinite_priority) : group->find_min();
    }

    Priority Delta = 0;
    Element no_element;
    Value no_value;
    Priority infinite_priority;
    PriorityQueue1<Element, Priority, Value> group_minima;
    std::vector<typename Group::ElementlementQueue::handle_type> elements;
};

} /* namespace Koala */
