#include <graph/GraphTools.hpp>

#include "cograph_rec/Cotree.hpp"

namespace Koala {
    Cotree::Cotree(NetworKit::Graph &Graph) {
        graph = &Graph;
        prepared = false;
    }

    void Cotree::BuildTree() {
        reverse(order.begin(), order.end());
        long long n = graph->numberOfNodes(), i;

        for (i = 0; i <= 2 * n; i++) {
            left_son.push_back(0);
            right_son.push_back(0);
            parent.push_back(0);
            type.push_back(0);
            size.push_back(0);
        }

        left_son[n] = order[0].first.first;
        right_son[n] = -1;
        parent[n] = -1;
        type[n] = 0;

        left_son[order[0].first.first] = -1;
        right_son[order[0].first.first] = -1;
        type[order[0].first.first] = 2;
        parent[order[0].first.first] = n;

        for (i = 1; i < n; i++) {
            left_son[order[i].first.first] = -1;
            right_son[order[i].first.first] = -1;
            type[order[i].first.first] = 2;
            parent[order[i].first.first] = n + i;

            long long p = parent[order[i].first.second];
            if (left_son[p] == order[i].first.second) {
                left_son[p] = n + i;
            } else {
                right_son[p] = n + i;
            }

            parent[order[i].first.second] = n + i;

            left_son[n + i] = order[i].first.first;
            right_son[n + i] = order[i].first.second;
            parent[n + i] = p;
            type[n + i] = order[i].second;
        }
        prepared = true;
    }
}