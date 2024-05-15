#include "pathwidth/Pathwidth.hpp"

namespace Koala {
    NetworKit::count Pathwidth::SubtreeSize(NetworKit::count n, NetworKit::count v) {
        if (recognition->cotree->nodes[v].type == 2) {
            recognition->cotree->nodes[v].size = 1;
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (recognition->cotree->nodes[v].left_son != -1) {
                l = SubtreeSize(n, recognition->cotree->nodes[v].left_son);
            }

            if (recognition->cotree->nodes[v].right_son != -1) {
                r = SubtreeSize(n, recognition->cotree->nodes[v].right_son);
            }
            recognition->cotree->nodes[v].size = r + l + 1;
            return r + l + 1;
        }
    }

    NetworKit::count Pathwidth::pathwidth(NetworKit::count n, NetworKit::count v) {
        if (recognition->cotree->nodes[v].type == 2) {
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (recognition->cotree->nodes[v].left_son != -1) {
                l = pathwidth(n, recognition->cotree->nodes[v].left_son);
            }

            if (recognition->cotree->nodes[v].right_son != -1) {
                r = pathwidth(n, recognition->cotree->nodes[v].right_son);
            }

            if (recognition->cotree->nodes[v].type == 0) {
                return std::max(l, r);
            } else {
                NetworKit::count minpathwidth = 4 * n;
                if (recognition->cotree->nodes[v].left_son != -1) {
                    minpathwidth = std::min(minpathwidth, r +
                                                          recognition->cotree->nodes[recognition->cotree->nodes[v].left_son].size);
                }

                if (recognition->cotree->nodes[v].right_son != -1) {
                    minpathwidth = std::min(minpathwidth, l +
                                                          recognition->cotree->nodes[recognition->cotree->nodes[v].right_son].size);
                }
                return minpathwidth;
            }
        }
    }

    void Pathwidth::run() {
        NetworKit::count n = graph->numberOfNodes();
        recognition = new Koala::CographRecognition(*graph);
        recognition->run();
        if (recognition->cotree->prepared == false) {
            recognition->cotree->BuildTree();
        }
        SubtreeSize(n, n);
        width = pathwidth(n, n);
    }
}