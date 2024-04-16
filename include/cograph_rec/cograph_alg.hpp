
#pragma once

#include <optional>
#include <vector>
#include <list>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace Koala {


    class CographRecognition : public NetworKit::Algorithm {
    public:
        enum class State {
            UNKNOWN,
            COGRAPH,
            IS_NOT_COGRAPH
        };
        State status;
        NetworKit::Graph *graph, *original_graph;

        long long num_of_parts, num_of_nodes;
        int prepared;

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

        part *l, *r, *first_part, *last_part;
        element *origin;
        std::vector<long long> left_son, right_son, parent, type, size;
        std::vector<std::pair<std::pair<long long, long long>, long long> > order;
        std::vector<element *> E;
        std::list<part *> unused_parts, garbage;


        bool IsCograph() const;

        CographRecognition::State GetState() const;

        explicit CographRecognition(NetworKit::Graph &graph);

        void run();

        void Check() const;

        void AddPart(part *l, part *r);

        void EraseElement(element *v);

        void ErasePart(part *p);

        void AddElementToPart(part *P, element *v);

        void Rcheck();

        void Lcheck();

        void Clear();

        void ShowTheOrder();

        void BuildTree();

        long long pathwidth(long long n, long long v);

        long long SubtreeSize(long long n, long long v);

        long long GetCographPathwidth();

        long long MaxCliqueSize(long long n, long long v);

        long long GetMaxCliqueSize();

        long long BruteForceCliqueSize();

    };

} /* namespace Koala */
