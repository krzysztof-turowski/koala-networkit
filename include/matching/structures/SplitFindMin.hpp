#pragma once

#include <list>
#include <vector>

namespace Koala {

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
template<typename Element, typename Cost, typename Value, typename Id>
class SplitFindMin {
 public:
    class List;
    /**
     * Creates a split-findmin data structure with level alpha(m, n)
     *
     * @param size maximum value of a stored element
     * @param infinite_cost default value for a cost
     * @param empty_val default value for a value associated with a cost
     * @param n parameter used in calculating the level
     * @param m parameter used in calculating the level
    */
    SplitFindMin(Element size, Cost infinite_cost, Value empty_val, int n, int m):
        SplitFindMin(size, infinite_cost, empty_val, alpha(m, n)) {}

    /**
     * Creates a list from a provided lists of elements with a given id.
     * All elements have default cost values.
     *
     * @param nodes elements to create a list of
     * @param id id associated with the list
     * @returns pointer to the new list
    */
    List* init(const std::list<Element>& nodes, Id id) {
        List* L = new List;
        L->id = id;
        L->i = max_i;
        L->sublist = nullptr;
        L->min_cost = infinite_cost;
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
    Id list(Element element) {
        return find_list(element, max_i)->id;
    }

    /**
     * Splits the list containing the provided element into two list - one containing all the
     * elements up to the provided one and a second one containg the remaining elements.
     * The two lists are given new ids.
     *
     * @param element the element on whose position to split it's list
     * @param id1 the new id of the list containg elements up to the provided one
     * @param id2 the new id of the list containg remaining elements
     * @returns pair of two lists resulting from the split
    */
    std::pair<List*, List*> split(Element element, Id id1, Id id2) {
        return split(element, id1, id2, max_i);
    }

    /**
     * Updates the cost of the element if the new one is lower
     *
     * @param element the element whose cost is to be updated
     * @param cost value of the cost
     * @param value value associated with the cost
    */
    void decreaseCost(Element element, Cost cost, Value value) {
        decrease_cost(element, cost, value, max_i);
    }

    /**
     * Returns the minimum cost of an element of the list and it's associated value
     *
     * @param L list whose minimum cost is returned
     * @returns pair containing the minimum cost in the list and it's associated value
    */
    std::pair<Cost, Value> findMin(List* L) {
        return {L->min_cost, L->min_val};
    }

    /**
     * Returns the current cost of an element and it's associated value
     *
     * @param element elements for which current cost and it's value is returned
     * @returns pair containing the current cost of the provided element and it's associated value
    */
    std::pair<Cost, Value> currentCost(Element element) {
        return {cost[max_i][element], val[max_i][element]};
    }

    /**
     * Deletes the list
     *
     * @param L the list to be deleted
    */
    void deleteList(List* L) {
        reset_costs(L);
        delete L;
    }

    struct Sublist {
     private:
        int level;
        bool head;
        List* list;
        List* elements;
        std::list<Sublist*>::iterator sublist_it;

        ~Sublist() { delete elements; }

        friend class List;
        friend class SplitFindMin;
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

        ~List() {
            clear();
        }

     private:
        Id id;
        int i;
        Cost min_cost;
        Value min_val;
        std::list<Element> nodes;
        std::list<Sublist*> head, tail;
        std::list<Element> head_singletons, tail_singletons;
        Sublist* sublist;

        friend class SplitFindMin;
    };

 private:
    // Hardcoded values of AcCostermann's function below 1000000000
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
        while (2 * A(i, j + 1) <= n) j++;
        return j;
    }

    static int alpha(int m, int n) {
        if (m < n) return 1;
        int i = 1;
        while (A(i, m / n) < n) i++;
        return i;
    }

    SplitFindMin(Element size, Cost infinite_cost, Value empty_val, int level):
        max_i(level),
        size(size),
        no_element(size + 1),
        infinite_cost(infinite_cost),
        empty_val(empty_val),
        e(level + 1),
        cost(level + 1),
        val(level + 1),
        element_list(level + 1),
        list_it(level + 1),
        superelement_nodes(level + 1) {
            for (int i = 0; i <= level; ++i) {
                e[i] = std::vector<Element>(size, no_element);
                element_list[i] = std::vector<List*>(size, nullptr);
                cost[i] = std::vector<Cost>(size, infinite_cost),
                val[i] = std::vector<Value>(size, empty_val),
                list_it[i] = std::vector<typename std::list<Element>::iterator>(size);
                superelement_nodes[i] = std::vector<std::list<Element>>(size);
            }
        }

    int max_i;
    Element size;
    Element no_element;
    Cost infinite_cost;
    Value empty_val;
    std::vector<std::vector<Element>> e;
    std::vector<std::vector<Cost>> cost;
    std::vector<std::vector<Value>> val;
    std::vector<std::vector<List*>> element_list;
    std::vector<std::vector<typename std::list<Element>::iterator>> list_it;
    std::vector<std::vector<std::list<Element>>> superelement_nodes;

    void initialize_head(List* L, const std::list<Element>& nodes, Id id, int i) {
        auto nodes_cpy = nodes;
        for (auto it = nodes_cpy.begin(); it != nodes_cpy.end(); ++it)
            list_it[i][*it] = it;
        L->nodes.splice(L->nodes.begin(), std::move(nodes_cpy));

        int j = a(i, nodes.size()) + 1;
        int remaining = nodes.size();
        auto it = nodes.rbegin();
        Sublist* sublist = nullptr;
        std::list<Element> sublist_elements;

        do {
            int prev_j = j;
            while (j >= 0 && 2 * A(i, j) > remaining) --j;

            if (j != prev_j) {
                if (sublist != nullptr) {
                    sublist->elements = new List;
                    sublist->elements->sublist = sublist;
                    sublist->elements->i = i - 1;
                    sublist->elements->id = id;
                    sublist->elements->min_cost = infinite_cost;
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
                Element singleton = *it;
                e[i][singleton] = no_element;
                element_list[i][singleton] = L;

                L->head_singletons.push_front(singleton);
                it++;
            } else {
                int size = 2 * A(i, j);
                Element superelement = *it;
                sublist_elements.push_front(superelement);
                cost[i-1][superelement] = infinite_cost;
                val[i-1][superelement] = empty_val;

                std::list<Element> se_nodes;
                for (int j = 0; j < size; ++j) {
                    se_nodes.push_front(*it);
                    it++;
                }

                for (auto n : se_nodes) {
                    e[i][n] = superelement;
                    update_cost_at_level(superelement, i - 1, cost[i][n], val[i][n]);
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
            sublist->elements->min_cost = infinite_cost;
            sublist->elements->min_val = empty_val;

            initialize_head(sublist->elements, sublist_elements, id, i - 1);

            L->head.push_front(sublist);
            sublist->sublist_it = L->head.begin();
        }

        for (auto n : nodes) update_cost(L, cost[i][n], val[i][n]);
    }

    void initialize_tail(List* L, const std::list<Element>& nodes, Id id, int i) {
        auto nodes_cpy = nodes;
        for (auto it = nodes_cpy.begin(); it != nodes_cpy.end(); ++it)
            list_it[i][*it] = it;
        L->nodes.splice(L->nodes.end(), std::move(nodes_cpy));

        int j = a(i, nodes.size()) + 1;
        int remaining = nodes.size();
        auto it = nodes.begin();
        Sublist* sublist = nullptr;
        std::list<Element> sublist_elements;

        do {
            int prev_j = j;
            while (j >= 0 && 2 * A(i, j) > remaining) --j;

            if (j != prev_j) {
                if (sublist != nullptr) {
                    sublist->elements = new List;
                    sublist->elements->sublist = sublist;
                    sublist->elements->id = id;
                    sublist->elements->i = i - 1;
                    sublist->elements->min_cost = infinite_cost;
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
                Element singleton = *it;
                e[i][singleton] = no_element;
                element_list[i][singleton] = L;
                L->tail_singletons.push_back(singleton);
                it++;
            } else {
                int size = 2 * A(i, j);

                std::list<Element> se_nodes;
                for (int j = 0; j < size; ++j) {
                    se_nodes.push_back(*it);
                    it++;
                }

                Element superelement = se_nodes.back();
                sublist_elements.push_back(superelement);
                cost[i-1][superelement] = infinite_cost;
                val[i-1][superelement] = empty_val;

                for (auto n : se_nodes) {
                    e[i][n] = superelement;
                    update_cost_at_level(superelement, i - 1, cost[i][n], val[i][n]);
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
            sublist->elements->min_cost = infinite_cost;
            sublist->elements->min_val = empty_val;

            initialize_tail(sublist->elements, sublist_elements, id, i - 1);

            L->tail.push_back(sublist);
            sublist->sublist_it = std::prev(L->tail.end());
        }

        for (auto n : nodes) update_cost(L, cost[i][n], val[i][n]);
    }

    List* find_list(Element u, int i) {
        if (is_singleton(u, i)) {
            return element_list[i][u];
        } else {
            auto L = find_list(e[i][u], i - 1);
            return L->sublist->list;
        }
    }

    List* decrease_cost(Element u, Cost x, Value v, int i) {
        update_cost(cost[i][u], val[i][u], x, v);
        if (is_singleton(u, i)) {
            auto L = element_list[i][u];
            update_cost(L, x, v);
            return L;
        } else {
            auto eL = decrease_cost(e[i][u], x, v, i - 1);
            auto L = eL->sublist->list;
            update_cost(L, x, v);
            return L;
        }
    }

    std::pair<List*, List*> split(Element x, Id id1, Id id2, int i) {
        if (is_singleton(x, i)) {
            List* L1 = element_list[i][x];
            L1->id = id1;
            List* L2 = new List;
            L2->id = id2;
            L2->i = i;
            L2->min_cost = infinite_cost;
            L2->min_val = empty_val;
            L2->nodes.splice(L2->nodes.begin(), L1->nodes,
                std::next(list_it[i][x]), L1->nodes.end());

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
            // Find superelemnt e(x)
            Element superelement = e[i][x];
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

            // S1 contains superlement before e(x) in the sublist
            // S2 contains superlement after e(x) in the sublist

            // Create a new sublist for S2
            Sublist* S2sublist = new Sublist;
            S2sublist->elements = S2;
            S2->sublist = S2sublist;
            S2sublist->level = S1->sublist->level;
            S2sublist->head = S1->sublist->head;

            // Split elements in superelment e(x) into those until and after x
            auto xit = std::find(superelement_nodes[i][superelement].begin(),
                                 superelement_nodes[i][superelement].end(), x);
            std::list<Element> nodes1(superelement_nodes[i][superelement].begin(), std::next(xit));
            std::list<Element> nodes2(std::next(xit), superelement_nodes[i][superelement].end());

            // Find the list
            auto L1 = S1->sublist->list;
            L1->id = id1;
            assert(L1->i == i);

            // Create a new list
            auto L2 = new List;
            L2->id = id2;
            L2->i = i;
            L2->min_cost = infinite_cost;
            L2->min_val = empty_val;

            // Update nodes lists for L1 and L2
            L2->nodes.splice(L2->nodes.begin(), L1->nodes,
                             std::next(list_it[i][superelement_nodes[i][superelement].back()]),
                             L1->nodes.end());

            L1->nodes.erase(list_it[i][superelement_nodes[i][superelement].front()],
                            L1->nodes.end());

            // Split the sublist list in L1
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

            // Update the cost of lists
            calculate_min(L1, i);
            calculate_min(L2, i);

            // Initialize the tail of L1 with nodes until x in e(x)
            initialize_tail(L1, nodes1, id1, i);

            // Initialize the head of L2 with nodes after x in e(x)
            if (nodes2.size() > 0)
                initialize_head(L2, nodes2, id2, i);

            return {L1, L2};
        }
    }

    void calculate_min(List* L, int i) {
        L->min_cost = infinite_cost;
        L->min_val = empty_val;

        for (auto hs : L->head_singletons) update_cost(L, cost[i][hs], val[i][hs]);
        for (auto sl : L->head) update_cost(L, sl->elements->min_cost, sl->elements->min_val);
        for (auto sl : L->tail) update_cost(L, sl->elements->min_cost, sl->elements->min_val);
        for (auto ts : L->tail_singletons) update_cost(L, cost[i][ts], val[i][ts]);
    }

    void set_list_pointers(List* L, int i) {
        for (auto hs : L->head_singletons) element_list[i][hs] = L;
        for (auto sl : L->head) sl->list = L;
        for (auto sl : L->tail) sl->list = L;
        for (auto ts : L->tail_singletons) element_list[i][ts] = L;
    }

    bool is_singleton(Element u, int i) {
        return e[i][u] == no_element;
    }

    void update_cost_at_level(Element u, int i, Cost c, Value v) {
        update_cost(cost[i][u], val[i][u], c, v);
    }

    void update_cost(List* L, Cost cost, Value val) {
        update_cost(L->min_cost, L->min_val, cost, val);
    }

    void update_cost(Cost& old_cost, Value& old_val, Cost cost, Value val) {
        if (cost < old_cost) {
            old_cost = cost;
            old_val = val;
        }
    }

    void reset_costs(List* L) {
        for (auto n : L->nodes) {
            cost[L->i][n] = infinite_cost;
            val[L->i][n] = empty_val;
        }
        for (auto s : L->head) reset_costs(s->elements);
        for (auto s : L->tail) reset_costs(s->elements);
    }
};

} /* namespace Koala */
