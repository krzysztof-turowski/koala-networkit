#pragma once

#include <graph/GraphTools.hpp>

namespace Koala {
    class Cotree {
    public:
        NetworKit::Graph *graph;
        std::vector<long long> left_son, right_son, parent, type, size;
        std::vector<std::pair<std::pair<long long, long long>, long long> > order;

        Cotree(NetworKit::Graph &Graph);

        bool prepared;

        void BuildTree();

        void SetOrder(std::vector<std::pair<std::pair<long long, long long>, long long> > a) {
            order = a;
        }

        long long pathwidth(long long n, long long v);

        long long SubtreeSize(long long n, long long v);

        long long GetCographPathwidth();

        long long MaxCliqueSize(long long n, long long v);

        long long GetMaxCliqueSize();

        long long BruteForceCliqueSize();

        long long MaxIndependetSetSize(long long n, long long v);

        long long GetMaxIndependetSetSize();

        long long BruteForceIndependetSetSize();
    };
}