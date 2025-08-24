#include <flow/minimum_cost_flow/OrlinMCF.hpp>


namespace Koala {

NetworKit::Graph::NodeIntAttribute OrlinMCF::base(NetworKit::Graph &g) {
    return g.getNodeIntAttribute("base");
}
    
NetworKit::Graph::NodeIntAttribute OrlinMCF::excess(NetworKit::Graph &g) {
    return g.getNodeIntAttribute("excess");
}

NetworKit::Graph::NodeIntAttribute OrlinMCF::potential(NetworKit::Graph &g) {
    return g.getNodeIntAttribute("potential");
}

NetworKit::Graph::EdgeIntAttribute OrlinMCF::flow(NetworKit::Graph &g) {
    return g.getEdgeIntAttribute("flow");
}

void OrlinMCF::runImpl() {

}

void OrlinMCF::initialize() {
    graph.attachNodeIntAttribute("base");
    graph.attachNodeIntAttribute("excess");
    graph.attachNodeIntAttribute("potential");

    graph.indexEdges();
    graph.attachEdgeIntAttribute("flow");

}

void OrlinMCF::makeConnected() {

}



void OrlinMCF::runImpl() {
    setupValues();
    NetworKit::Graph contGraph = graph;
    auto flowAttr = flow(contGraph);
    auto excessAttr = excess(contGraph);


    while (imbalanced) {
        bool zero = true; 
        for (auto f : flowAttr) {
            if (f.second) {
                zero = false;
            } 
        }
        int maxDelta = 0;
        if (zero) { 
            for (auto ex : excessAttr) {
                maxDelta = std::max(ex.second, maxDelta);
            }
        }
        delta = std::min(delta, maxDelta);

        deltaScalingPhase(contGraph);
        delta /= 2;
    }
}


void OrlinMCF::deltaScalingPhase(NetworKit::Graph &g) {
    contractionPhase(g);
    auto excessAttr = excess(g);

    std::vector<NetworKit::node> s,t;
    g.forNodes([&](NetworKit::node node) {
        int ex = excessAttr.get(node);
        if (ex >= ALPHA*delta) {
            s.push_back(node);
        }
        if (ex <= -ALPHA*delta) {
            t.push_back(node);
        }
    });

    while (s.size() && t.size()) {
        NetworKit::node k = s.back();
        NetworKit::node v = t.back();

        // compute d
        

    
    }
        
}

void OrlinMCF::contractionPhase(NetworKit::Graph &g) {
    auto pi = potential(g);

    auto forAllNeighbours = [&](std::function<void(const NetworKit::node, const int)> const& fun) {
        for (auto v: g.nodeRange()) {
            for (int nei = 0; nei < g.degreeOut(v); ++nei) {
                fun(v, nei);
            }
        }
    };
    
    auto isBigArc = [&]() -> 
        std::tuple<bool, NetworKit::node, NetworKit::node> {
        auto f = flow(g);
        for (auto v: g.nodeRange()) {
            for (int nei = 0; nei < g.degreeOut(v); ++nei) {
                auto [w, id] = g.getIthNeighborWithId(v, nei);
                if (f.get(id) > 3*graph.numberOfNodes()*delta) {
                    return { true, v,w };
                }
            }
        }
        return { false, 0, 0 };
    };
    
    while (true) {
        auto [is, v, w] = isBigArc();
        if (!is) {
            break;
        }
        
        contractNodes(g, v, w);
        
        forAllNeighbours(
            [&](NetworKit::node v, int neigh) {
                auto [w, weight] = g.getIthNeighborWithWeight(v, neigh);
                int pi_v = pi.get(v);
                int pi_w = pi.get(w);
                g.setWeightAtIthInNeighbor(NetworKit::unsafe, v, neigh, weight - pi_v + pi_w);
            }
        );
        for(auto p : pi) {
            pi.set(p.first, 0);
        }
    }


}

void OrlinMCF::setupValues() {
    auto baseAttr = base(graph);
    auto excessAttr = excess(graph);
    auto potentialAttr = potential(graph);

    imbalanced = 0;
    delta = 0;

    graph.forNodes(
        [&](NetworKit::node v) {
            int e = baseAttr.get(v);
            excessAttr.set(v, e);
            potentialAttr.set(v, 0);
            if (e) {
                ++imbalanced;
            }

            delta = std::max(e, delta);
        }
    );
}

void OrlinMCF::contractNodes(NetworKit::Graph &g, NetworKit::node v, NetworKit::node w) {
    auto flowAttr = flow(g);
    if (w < v) {
        std::swap(v, w);
    }

    g.forEdgesOf(w, 
        [&](NetworKit::node, NetworKit::node y, NetworKit::edgeweight weight, NetworKit::edgeid id) {
            if(y == v) return;
            graph.addEdge(v, y, weight);
            int upper = graph.upperEdgeIdBound();
            flowAttr.set(upper - 1, flowAttr.get(id));
    });
    g.forInEdgesOf(w,
        [&](NetworKit::node, NetworKit::node y, NetworKit::edgeweight weight, NetworKit::edgeid id) {
            if(y == v) 
                return;
            graph.addEdge(y, v, weight);
            int upper = graph.upperEdgeIdBound();
            flowAttr.set(upper - 1, flowAttr.get(id));
    });
    
    
    auto ex = excess(g);
    auto bs = base(g);

    ex.set(v, ex.get(w) + ex.get(v));
    bs.set(v, bs.get(w) + bs.get(v));

    g.removeNode(w);
}

} /* namespace Koala */