#include <matching/MaximumMatching.hpp>

namespace Koala {

EdmondsMaximumMatching::EdmondsMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        current_blossom(graph.upperNodeIdBound(), nullptr),
        y(graph.upperNodeIdBound(), max_weight / 2) { 

    graph.forNodes([this] (NetworKit::node vertex) {
        current_blossom[vertex] = trivial_blossom[vertex];
    });
}

void EdmondsMaximumMatching::scan_edges(Blossom* b) {
    for (auto v : b->nodes) {
        graph.forEdgesOf(v, [this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            if (is_tight(u, v, id)) 
                edge_queue.push({u, v, id});
        });
    };
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
    edge_queue = {};

    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        if (is_tight(u, v, id)) {
            edge_queue.push({u, v, id});
        } else if (is_tight(v, u, id)) {
            edge_queue.push({v, u, id});
        }
    });
}

bool EdmondsMaximumMatching::has_useful_edges() {
    return !edge_queue.empty();
}

EdmondsMaximumMatching::Edge EdmondsMaximumMatching::get_useful_edge() {
    auto edge = edge_queue.front(); edge_queue.pop();
    return edge;
}

void EdmondsMaximumMatching::handle_grow(Blossom* odd_blossom, Blossom* even_blossom) {
    scan_edges(even_blossom);
}

void EdmondsMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    for (auto [b, edge] : new_blossom->subblossoms) { 
        for (auto v : b->nodes)
            current_blossom[v] = new_blossom;

        if (b->label == odd)
            scan_edges(b);
    }
}

void EdmondsMaximumMatching::handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) {}

void EdmondsMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        for (auto v : b->nodes)
            current_blossom[v] = b;

        if (b->label == even)
            scan_edges(b);
    }
}

void EdmondsMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        for (auto v : b->nodes)
            current_blossom[v] = b;
    }
}

void EdmondsMaximumMatching::adjust_by_delta(MaximumWeightMatching::weight delta) {
    graph.forNodes([this, delta] (NetworKit::node v) {
        Blossom* v_blossom = get_blossom(v);
        if (v_blossom->label == even) {
            y[v] -= delta;
        } else if (v_blossom->label == odd) {
            y[v] += delta;
        }
    });

    for (Blossom* blossom : blossoms) {
        if (blossom->label == even) {
            blossom->z += 2 * delta;
        } else if (blossom->label == odd) {
            blossom->z -= 2 * delta;
        }
    }
}

MaximumWeightMatching::weight EdmondsMaximumMatching::calc_delta1() {
    // min u_i : i - even vertex

    MaximumWeightMatching::weight res = infinite_weight;
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (get_blossom(v)->label == even) {
            res = std::min(res, y[v]);
        }
    });
    return res;
}

MaximumWeightMatching::weight EdmondsMaximumMatching::calc_delta2() {
    // min pi_ij : i - even vertex, j - free vertex

    MaximumWeightMatching::weight res = infinite_weight;
    graph.forEdges([this, &res] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        Blossom* u_blossom = get_blossom(u);
        Blossom* v_blossom = get_blossom(v);
        if ((u_blossom->label == even && v_blossom->label == free) || 
            (u_blossom->label == free && v_blossom->label == even)) {
            res = std::min(res, slack(id));
        }
    });
    return res;
}

MaximumWeightMatching::weight EdmondsMaximumMatching::calc_delta3() {
    // min pi_ij / 2 : i,j - even vertices in different blossoms

    MaximumWeightMatching::weight res = infinite_weight;
    graph.forEdges([this, &res] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        Blossom* u_blossom = get_blossom(u);
        Blossom* v_blossom = get_blossom(v);
        if (u_blossom != v_blossom && u_blossom->label == even && v_blossom->label == even) {
            res = std::min(res, slack(id) / 2);
        }
    });
    return res;
}

MaximumWeightMatching::weight EdmondsMaximumMatching::calc_delta4() {
    // min z_k / 2 : B_k - odd blossom 

    MaximumWeightMatching::weight res = infinite_weight;
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

MaximumWeightMatching::weight EdmondsMaximumMatching::slack(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    Blossom* u_blossom = get_blossom(u);
    Blossom* v_blossom = get_blossom(v);
    return y[u] + y[v] - w + (u_blossom == v_blossom ? u_blossom->z : 0);
}

bool EdmondsMaximumMatching::is_tight(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge) {
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);
    return !is_in_matching[edge] && slack(edge) == 0 && 
            u_blossom->label == even && 
            ((v_blossom->label == even && u_blossom != v_blossom) || v_blossom->label == free);
}

void EdmondsMaximumMatching::check_consistency() {
    std::cerr << "Current vertices: \n";
    graph.forNodes([this] (NetworKit::node v) {
        if (this->matched_vertex[v] != NetworKit::none)
            std::cerr << v << " " << this->matched_vertex[v] << " : ";
        else std::cerr << v << " - : ";
        get_blossom(v)->nodes_print(); 
        std::cerr << std::endl;
    });
    std::cerr << "Current weights: \n";
    for (int i = 0; i < y.size(); ++ i) {
        std::cerr << i << ": " << y[i] << std::endl;
    }
    for (auto b : blossoms) {
        if (!b->is_trivial()) {
            b->short_print(); 
            std::cerr << " : " << b->z << std::endl; 
        }
    }
    std::cerr << "Edge queue:\n";
    std::queue<Edge> q = edge_queue;
    while (!q.empty()) {
        std::cerr << edge_to_string(q.front()) << std::endl;
        q.pop();
    }
    std::cerr << "Edges:\n";
    for (NetworKit::edgeid id = 0; id < graph.upperEdgeIdBound(); ++ id) {
        auto [u, v, w] = graph_edges[id];
        std::cerr << edge_to_string({u, v, id}) << " : " << slack(id) << std::endl; 
    }
}

} /* namespace Koala */
