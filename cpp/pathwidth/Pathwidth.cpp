#include "pathwidth/Pathwidth.hpp"

namespace Koala {
    NetworKit::count Pathwidth::subtree_size(NetworKit::count v) {
        CoNode V(0,0,0);
        V=cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            V.size = 1;
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (V.left_son != NetworKit::none) {
                l = subtree_size(V.left_son);
            }

            if (V.right_son != NetworKit::none) {
                r = subtree_size(V.right_son);
            }
            V.size = r + l + 1;
            return r + l + 1;
        }
    }

    NetworKit::count Pathwidth::pathwidth(NetworKit::count v) {
        CoNode V(0,0,0);
        V=cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (V.left_son != NetworKit::none) {
                l = pathwidth(V.left_son);
            }

            if (V.right_son != NetworKit::none) {
                r = pathwidth(V.right_son);
            }

            if (V.type == NodeType::UNION_NODE) {
                return std::max(l, r);
            } else {
                NetworKit::count minpathwidth = 1000000 ;
                if (V.left_son != NetworKit::none) {
                    minpathwidth = std::min(minpathwidth, r +cotree.getNode(V.left_son).size);
                }

                if (V.right_son != NetworKit::none) {
                    minpathwidth = std::min(minpathwidth, l +cotree.getNode(V.right_son).size);
                }
                return minpathwidth;
            }
        }
    }

    void Pathwidth::run() {
        if (cotree.prepared == false) {
            cotree.buildTree();
        }
        subtree_size(cotree.graph->numberOfNodes());
        width = pathwidth(cotree.graph->numberOfNodes());
    }
}