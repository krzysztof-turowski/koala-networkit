//
// Created by scales on 16.04.24.
//
#include <list>
#include <set>

#include <unordered_map>
#include <map>

#include <graph/GraphTools.hpp>

std::vector<int> used;


void StartTest(long long n) {
    for (int i; i < n; i++) {
        used.push_back(0);
    }
}

bool FalseTwins(long long A, long long B, NetworKit::Graph *graph, long long twins_counter) {
    twins_counter += 2;
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


bool TrueTwins(long long A, long long B, NetworKit::Graph *graph, long long twins_counter) {
    twins_counter += 2;
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

int Twins(long long A, long long B, NetworKit::Graph *graph, long long twins_counter) {
    if (A == -1 || B == -1) {
        return 2;
    }
    if (TrueTwins(A, B, graph, twins_counter)) {
        return 1;
    }
    if (FalseTwins(A, B, graph, twins_counter)) {
        return 0;
    }
    return 2;
}