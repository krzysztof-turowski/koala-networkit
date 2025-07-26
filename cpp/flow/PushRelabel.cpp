#include <limits>
#include <algorithm>
#include <flow/PushRelabel.hpp>

using edge = NetworKit::Edge;

namespace Koala {

edge PushRelabel::reverse(const edge &p) {
    return NetworKit::Edge(p.v, p.u);
}
void PushRelabel::initialize() {
    V = graph->numberOfNodes();
    graph->forNodes([&](NetworKit::node v) {
        height[v] = 0;
        nextedge[v] = 0;
        excess[v] = 0;
    });
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        auto p = NetworKit::Edge(u, v);
        if (graph->addEdge(v, u, 0, true)) {
            capacity[reverse(p)] = 0;
        }
        flow[p] = 0;
        flow[reverse(p)] = 0;
        capacity[p] = w;
    });
    height[source] = V;
    excess[source] = std::numeric_limits<int>::max();
}

void PushRelabel::push(const NetworKit::node& vertex, const edge& e) {
    int fl = std::min(capacity[e] - flow[e], excess[vertex]);
    excess[vertex] -= fl;
    excess[e.v] += fl;
    flow[e] += fl;
    flow[reverse(e)] -= fl;

    if (excess[e.v] > 0 && e.v != source && e.v != target) {
        q.push(e.v);
    }
}

void PushRelabel::relabel(const NetworKit::node& vertex) {
    int minimum = std::numeric_limits<int>::max();

    graph->forNeighborsOf(vertex, [&](NetworKit::node w) {
        auto e = NetworKit::Edge(vertex, w);
        if (capacity[e] - flow[e] > 0) {
            minimum = std::min(minimum, height[e.v]);
        }
    });

    if (minimum != std::numeric_limits<int>::max()) {
        height[vertex] = minimum + 1;
    }
}

void PushRelabel::discharge(const NetworKit::node& vertex) {
    while (excess[vertex] > 0) {
        while (nextedge[vertex] < graph->degreeOut(vertex)) {
            auto u = graph->getIthNeighbor(vertex, nextedge[vertex]);
            auto e = NetworKit::Edge(vertex, u);

            if (capacity[e] - flow[e] > 0 && height[vertex] == height[e.v] + 1) {
                push(vertex, e);
                if (excess[vertex] == 0) {
                    return;
                }
            }
            nextedge[vertex]++;
        }

        relabel(vertex);
        nextedge[vertex] = 0;

        if (height[vertex] >= V) break;
    }
}

void PushRelabel::run() {
    initialize();

    graph->forNeighborsOf(source, [&](NetworKit::node w) {
        push(source, NetworKit::Edge(source, w));
    });

    while (!q.empty()) {
        auto vertex = q.front();
        q.pop();
        discharge(vertex);
    }

    flow_size = excess[target];
    hasRun = true;
}
}  // namespace Koala
