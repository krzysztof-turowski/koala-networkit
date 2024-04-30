#pragma once

#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {
    class element;

    class part {
    public:
        element *first, *last, *pivot;
        part *previous, *next, *division;
        long long size, amount;

        part() {
            size = 0;
            amount = 0;
            first = nullptr;
            last = nullptr;
            pivot = nullptr;
            previous = nullptr;
            next = nullptr;
            division = nullptr;
        }
    };

    class element {
    public:
        part *my_part;
        element *previous, *next;
        long long num = -1;
    };
};


