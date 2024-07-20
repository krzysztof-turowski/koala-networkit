#include <graph/GraphTools.hpp>

#include "recognition/сograph/Twins.hpp"

namespace Koala {

Twins::Twins(NetworKit::Graph &Graph) : graph(Graph) {
    used.resize(graph.numberOfNodes(), 0);
}

bool Twins::false_twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter) {
    NetworKit::count counter = 0;
    if (graph.degree(A) != graph.degree((B))) {
        return false;
    }

    for (auto node1 : graph.neighborRange(A)) {
        used[node1] = twins_counter;
        counter++;
    }

    for (auto node2 : graph.neighborRange(B)) {
        if (used[node2] == twins_counter) {
            counter--;
        } else {
            counter += graph.numberOfNodes();
        }
    }
    return (counter == 0);
}

bool Twins::true_twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter) {
    NetworKit::count counter = 0, flag = 0;
    if (graph.degree(A) != graph.degree((B))) {
        return false;
    }
    for (auto node1 : graph.neighborRange(A)) {
        if (node1 != B) {
            used[node1] = twins_counter;
            counter++;
        } else {
            flag++;
        }
    }
    for (auto node2 : graph.neighborRange(B)) {
        if (node2 != A) {
            if (used[node2] == twins_counter) {
                counter--;
            } else {
                counter += graph.numberOfNodes();
            }
        } else {
            flag++;
        }
    }
    return (counter == 0 && flag == 2);
}

int Twins::twins(NetworKit::count A, NetworKit::count B, NetworKit::count twins_counter) {
    if (A == NetworKit::none || B == NetworKit::none) {
        return 2;
    }
    twins_counter += 2;
    if (true_twins(A, B, twins_counter)) {
        return 1;
    }
    twins_counter += 2;
    if (false_twins(A, B, twins_counter)) {
        return 0;
    }
    return 2;
}

} /* namespace Koala */
