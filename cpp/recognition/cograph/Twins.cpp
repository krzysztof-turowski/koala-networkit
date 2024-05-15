#include <graph/GraphTools.hpp>

namespace Koala {
    std::vector<int> used;

    void StartTest(NetworKit::count n) {
        used.clear();
        for (int i = 0; i < n; i++) {
            used.push_back(0);
        }
    }

    bool FalseTwins(NetworKit::count A, NetworKit::count B, NetworKit::Graph *graph, NetworKit::count twins_counter) {
        int counter = 0;
        if (graph->degree(A) != graph->degree((B))) {
            return false;
        }

        for (int node1: graph->neighborRange(A)) {
            used[node1] = twins_counter;
            counter++;
        }

        for (auto node2: graph->neighborRange(B)) {
            if (used[node2] == twins_counter) {
                counter--;
            } else {
                counter += graph->numberOfNodes();
            }
        }
        return (counter == 0);
    }

    bool TrueTwins(NetworKit::count A, NetworKit::count B, NetworKit::Graph *graph, NetworKit::count twins_counter) {
        int counter = 0, flag = 0;
        if (graph->degree(A) != graph->degree((B))) {
            return false;
        }

        for (int node1: graph->neighborRange(A)) {
            if (node1 != B) {
                used[node1] = twins_counter;
                counter++;
            } else {
                flag++;
            }
        }

        for (auto node2: graph->neighborRange(B)) {
            if (node2 != A) {
                if (used[node2] == twins_counter) {
                    counter--;
                } else {
                    counter += graph->numberOfNodes();
                }
            } else {
                flag++;
            }
        }
        return (counter == 0 && flag == 2);
    }

    int Twins(NetworKit::count A, NetworKit::count B, NetworKit::Graph *graph, NetworKit::count twins_counter) {
        if (A == -1 || B == -1) {
            return 2;
        }
        twins_counter += 2;
        if (TrueTwins(A, B, graph, twins_counter)) {
            return 1;
        }
        twins_counter += 2;
        if (FalseTwins(A, B, graph, twins_counter)) {
            return 0;
        }
        return 2;
    }
}