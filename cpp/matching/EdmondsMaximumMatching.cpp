#include <matching/MaximumMatching.hpp>

namespace Koala {

EdmondsMaximumMatching::EdmondsMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        current_blossom(graph.upperNodeIdBound(), nullptr) { 
    
    MaximumMatching::edgeweight  max_weight = std::numeric_limits<MaximumMatching::edgeweight >::min();
    for (auto [u, v, w] : graph_edges)
        max_weight = std::max(w, max_weight);
    U = std::vector<MaximumMatching::edgeweight>(graph.upperNodeIdBound(), max_weight / 2);
    graph.forNodes([this] (NetworKit::node vertex) {
        current_blossom[vertex] = trivial_blossom[vertex];
    });
}

void EdmondsMaximumMatching::initialize_stage() {
    // Start search in all exposed blossoms 
    for (auto blossom : blossoms) {
        blossom->label = is_exposed(blossom) ? even : free;
        blossom->backtrack_edge = {NetworKit::none, NetworKit::none, NetworKit::none};
    }
}

void EdmondsMaximumMatching::finish_stage() {
    // Expand even blossoms with B_z = 0
    std::vector<Blossom*> to_expand;
    for (auto b : blossoms) {
        if (b->label == even && b->z == 0.0 && !b->is_trivial()) {
            to_expand.push_back(b);
        }
    }
    for (auto b : to_expand) expand_even_blossom(b);
}

void EdmondsMaximumMatching::initialize_substage() {
    useful_edges = {};

    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        if (is_useful(u, v, id)) {
            useful_edges.push({u, v, id});
        } else if (is_useful(v, u, id)) {
            useful_edges.push({v, u, id});
        }
    });
}

bool EdmondsMaximumMatching::has_useful_edges() {
    return !useful_edges.empty();
}

EdmondsMaximumMatching::EdgeInfo EdmondsMaximumMatching::get_useful_edge() {
    auto edge = useful_edges.front(); useful_edges.pop();
    return edge;
}

void EdmondsMaximumMatching::label_odd(Blossom* b) {}

void EdmondsMaximumMatching::label_even(Blossom* b) {
    // Check for new useful edges
    b->for_nodes([this] (NetworKit::node vertex) {
        graph.forEdgesOf(vertex, [this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            if (is_useful(u, v, id)) 
                useful_edges.push({u, v, id});
        });
    });
}

void EdmondsMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    for (auto [b, edge] : new_blossom->subblossoms) { 
        b->for_nodes([this, new_blossom] (NetworKit::node v) {
            current_blossom[v] = new_blossom;
        });
        if (b->label == odd) {
            label_even(b);
        }
    }
}

void EdmondsMaximumMatching::handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) {}

void EdmondsMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        Blossom* _b = b;
        b->for_nodes([this, _b] (NetworKit::node v) {
            current_blossom[v] = _b;
        });
        if (b->label == even) {
            label_even(b);
        }
    }
}

void EdmondsMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        Blossom* _b = b;
        b->for_nodes([this, _b] (NetworKit::node v) {
            current_blossom[v] = _b;
        });
    }
}

void EdmondsMaximumMatching::adjust_by_delta(MaximumMatching::edgeweight  delta) {
    graph.forNodes([this, delta] (NetworKit::node v) {
        Blossom* v_blossom = get_blossom(v);
        if (v_blossom->label == even) {
            U[v] -= delta;
        } else if (v_blossom->label == odd) {
            U[v] += delta;
        }
    });

    for (Blossom* blossom : blossoms) {
        if (blossom->label == even) {
            blossom->z += 2.0 * delta;
        } else if (blossom->label == odd) {
            blossom->z -= 2.0 * delta;
        }
    }
}

MaximumMatching::edgeweight  EdmondsMaximumMatching::calc_delta1() {
    // min u_i : i - even vertex

    MaximumMatching::edgeweight  res = std::numeric_limits<MaximumMatching::edgeweight >::max();
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (get_blossom(v)->label == even) {
            res = std::min(res, U[v]);
        }
    });
    return res;
}

MaximumMatching::edgeweight  EdmondsMaximumMatching::calc_delta2() {
    // min pi_ij : i - even vertex, j - free vertex

    MaximumMatching::edgeweight  res = std::numeric_limits<MaximumMatching::edgeweight >::max();
    graph.forEdges([this, &res] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        Blossom* u_blossom = get_blossom(u);
        Blossom* v_blossom = get_blossom(v);
        if ((u_blossom->label == even && v_blossom->label == free) || 
            (u_blossom->label == free && v_blossom->label == even)) {
            res = std::min(res, edge_dual_variable(id));
        }
    });
    return res;
}

MaximumMatching::edgeweight  EdmondsMaximumMatching::calc_delta3() {
    // min pi_ij / 2 : i,j - even vertices in different blossoms

    MaximumMatching::edgeweight  res = std::numeric_limits<MaximumMatching::edgeweight >::max();
    graph.forEdges([this, &res] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        Blossom* u_blossom = get_blossom(u);
        Blossom* v_blossom = get_blossom(v);
        if (u_blossom != v_blossom && u_blossom->label == even && v_blossom->label == even) {
            res = std::min(res, edge_dual_variable(id) / 2);
        }
    });
    return res;
}

MaximumMatching::edgeweight  EdmondsMaximumMatching::calc_delta4() {
    // min z_k / 2 : B_k - odd blossom 

    MaximumMatching::edgeweight  res = std::numeric_limits<MaximumMatching::edgeweight >::max();
    for (Blossom* b : blossoms) {
        if (b->label == odd && !b->is_trivial()) {
            res = std::min(res, b->z / 2);
        }
    }

    return res;
}

void EdmondsMaximumMatching::find_delta2_useful_edges() {}
void EdmondsMaximumMatching::find_delta3_useful_edges() {}

std::vector<EdmondsMaximumMatching::Blossom*> EdmondsMaximumMatching::get_odd_blossoms_to_expand() {
    std::vector<Blossom*> to_expand;
    for (auto b : blossoms) {
        if (b->label == odd && b->z == 0.0 && !b->is_trivial()) {
            to_expand.push_back(b);
        }
    }
    return to_expand;
}

EdmondsMaximumMatching::Blossom* EdmondsMaximumMatching::get_blossom(NetworKit::node vertex) {
    return current_blossom[vertex];
}

MaximumMatching::edgeweight  EdmondsMaximumMatching::edge_dual_variable(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    Blossom* u_blossom = get_blossom(u);
    Blossom* v_blossom = get_blossom(v);
    return U[u] + U[v] - w + (u_blossom == v_blossom ? u_blossom->z : 0);
}

bool EdmondsMaximumMatching::is_useful(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge) {
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);
    return !is_in_matching[edge] && edge_dual_variable(edge) == 0 && 
            u_blossom->label == even && 
            ((v_blossom->label == even && u_blossom != v_blossom) || v_blossom->label == free);
}

void EdmondsMaximumMatching::check_consistency() {
    std::cerr << "Current vertices: \n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << ": ";
        get_blossom(v)->short_print(); 
        std::cerr << std::endl;
    });
    std::cerr << "Current weights: \n";
    for (int i = 0; i < U.size(); ++ i) {
        std::cerr << i << ": " << U[i] << std::endl;
    }
    for (auto b : blossoms) {
        if (!b->is_trivial()) {
            std::cerr << b->z << " "; 
            b->short_print(); 
            std::cerr << std::endl; 
        }
    }
    std::cerr << "Edge queue:\n";
    std::queue<EdgeInfo> q = useful_edges;
    while (!q.empty()) {
        auto [u, v, id] = q.front(); q.pop();
        std::cerr << "(" << u << ", " << v << ")\n";
    }
    std::cerr << "Edges:\n";
    for (int id = 0; id < graph.upperEdgeIdBound(); ++ id) {
        auto [u, v, w] = graph_edges[id];
        std::cerr << "(" << u << ", " << v << ") : " << edge_dual_variable(id) << std::endl; 
    }
}

} /* namespace Koala */
