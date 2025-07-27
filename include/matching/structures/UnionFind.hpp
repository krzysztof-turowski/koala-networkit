#pragma once

#include <vector>

namespace Koala {

/**
 * Implementation of the union-find data structure operating on elements from a provided
 * interval. Allows for linking sets of two provided elements.
 * Each set is associated with a value that can be set when linking two sets and that is returned
 * when calling the find method.
 *
 * @tparam Element type of elements
 * @tparam Id type of the ID assigned to each set
*/
template <typename Element, typename Id>
class UnionFind {
 public:
    /**
     * Creates a union-find data structure for storing elements between 0 and size - 1
     *
     * @param size upper bound on the elements that can be stored
    */
    explicit UnionFind(Element size): size(size), root(size), parent(size), rank(size) {}

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
    void reset(Element element, Id id) {
        parent[element] = element;
        root[element] = id;
        rank[element] = 0;
    }

    /**
     * Links sets of two provided elements.
     *
     * @param x first element
     * @param y second element
     * @param value the value that is to associated with the created set
    */
    void link(Element x, Element y, Id value) {
        Element rx = find_root(x);
        Element ry = find_root(y);

        if (rank[rx] < rank[ry]) {
            parent[rx] = ry;
            root[ry] = value;
        } else {
            parent[ry] = rx;
            root[rx] = value;
            if (rank[rx] == rank[ry]) {
                rank[rx]++;
            }
        }
    }

    /**
     * Returns the ID of the set that the provided element belongs to.
     *
     * @param element value of the element
    */
    Id find(Element element) {
        return root[find_root(element)];
    }

 private:
    Element size;
    std::vector<Id> root;
    std::vector<Element> parent;
    std::vector<int> rank;

    Element find_root(Element x) {
        return parent[x] == x ? x : (parent[x] = find_root(parent[x]));
    }
};

} /* namespace Koala */
