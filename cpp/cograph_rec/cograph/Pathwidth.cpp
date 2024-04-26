
//
// Created by scales on 16.04.24.
//
#include <list>
#include <set>
#include <unordered_map>
#include <map>
#include <graph/GraphTools.hpp>

#include "cograph_rec/Cotree.hpp"

namespace Koala {

    long long Cotree::SubtreeSize(long long n, long long v) {
        if (type[v] == 2) {
            size[v] = 1;
            return 1;
        } else {
            long long l = 0, r = 0;
            if (left_son[v] != -1) {
                l = SubtreeSize(n, left_son[v]);
            }

            if (right_son[v] != -1) {
                r = SubtreeSize(n, right_son[v]);
            }
            size[v] = r + l + 1;
            return r + l + 1;
        }
    }

    long long Cotree::pathwidth(long long n, long long v) {
        if (type[v] == 2) {
            return 1;
        } else {
            long long l = 0, r = 0;

            if (left_son[v] != -1) {
                l = pathwidth(n, left_son[v]);
            }

            if (right_son[v] != -1) {
                r = pathwidth(n, right_son[v]);
            }

            if (type[v] == 0) {
                return std::max(l, r);
            } else {
                long long minpathwidth = 4 * n;
                if (left_son[v] != -1) {
                    minpathwidth = std::min(minpathwidth, r + size[left_son[v]]);
                }

                if (right_son[v] != -1) {
                    minpathwidth = std::min(minpathwidth, l + size[right_son[v]]);
                }
                return minpathwidth;
            }
        }
    }


    long long Cotree::GetCographPathwidth() {
        long long n = graph->numberOfNodes();
        if (prepared == false) {
            BuildTree();
        }
        SubtreeSize(n, n);
        return pathwidth(n, n);
    }
}