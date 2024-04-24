#include <matching/MaximumMatching.hpp>

namespace Koala {

GabowMaximumMatching::GabowMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        current_blossom(graph.upperNodeIdBound(), nullptr),
        U(graph.upperNodeIdBound(), 0),
        best_edge(graph.upperNodeIdBound(), no_edge) { 
    
    graph.sortEdges();
    MaximumMatching::edgeweight max_weight = std::numeric_limits<MaximumMatching::edgeweight>::min();
    graph.forEdges([&max_weight] 
            (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) { 
        max_weight = std::max(static_cast<MaximumMatching::edgeweight>(w), max_weight);
    });
    U = std::vector<MaximumMatching::edgeweight>(graph.upperNodeIdBound(), max_weight / 2);
    graph.forNodes([this] (NetworKit::node vertex) {
        this->current_blossom[vertex] = this->trivial_blossom[vertex];
    });

    for (auto b : trivial_blossom) {
        b->data = new GabowBlossomData();
    }
}

GabowMaximumMatching::GabowBlossomData* GabowMaximumMatching::get_data(Blossom* b) {
    return static_cast<GabowBlossomData*>(b->data);
}

void GabowMaximumMatching::initialize_stage() {
    edge_queue = {};
    // Start search in all exposed blossoms 
    graph.forNodes([this] (NetworKit::node v) {
        best_edge[v] = no_edge;
    });

    for (auto blossom : blossoms) {
        blossom->label = is_exposed(blossom) ? even : free;
        GabowBlossomData* data = get_data(blossom);
        data->best_edges = {};
        blossom->backtrack_edge = data->best_edge = no_edge;
    }
    for (auto blossom : blossoms) {
        if (blossom->label == even) {
            calc_best_edges(blossom);
        }
    }
}

void GabowMaximumMatching::finish_stage() {
    // Expand even blossoms with B_z = 0
    std::vector<Blossom*> to_expand;
    for (auto b : blossoms) {
        if (b->label == even && b->z == 0.0 && !b->is_trivial()) {
            to_expand.push_back(b);
        }
    }
    for (auto b : to_expand) expand_even_blossom(b);
}

void GabowMaximumMatching::initialize_substage() {}

bool GabowMaximumMatching::has_useful_edges() {
    return !edge_queue.empty();
}

GabowMaximumMatching::EdgeInfo GabowMaximumMatching::get_useful_edge() {
    auto edge = edge_queue.front(); edge_queue.pop();
    return edge;
}

void GabowMaximumMatching::label_odd(Blossom* b) {}

void GabowMaximumMatching::label_even(Blossom* b) {
    // Check for new useful edges
    calc_best_edges(b);
}

void GabowMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    for (auto [b, edge] : new_blossom->subblossoms) { 
        b->for_nodes([this, new_blossom] (NetworKit::node v) {
            this->current_blossom[v] = new_blossom;
        });
    }
    GabowBlossomData* data = new GabowBlossomData();
    new_blossom->data = data;
    auto best_slack = edge_slack(data->best_edge.id);
    for (auto [b, edge] : new_blossom->subblossoms) {
        calc_best_edges(b);
        auto b_data = get_data(b);
        for (auto [v_blossom, edge] : b_data->best_edges) {
            auto [u, v, id] = edge;
            if (get_blossom(v) == new_blossom) continue;
            auto slack = edge_slack(edge.id);
            if (data->best_edges.find(v_blossom) == data->best_edges.end() 
                    || slack < edge_slack(data->best_edges[v_blossom].id)) {
                data->best_edges[v_blossom] = edge;
                if (slack < best_slack) {
                    best_slack = slack;
                    data->best_edge = edge;
                }
            }
        }
        b_data->best_edges.clear();
        b_data->best_edge = no_edge;
    }
}

void GabowMaximumMatching::handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) {}

void GabowMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        Blossom* _b = b;
        b->for_nodes([this, _b] (NetworKit::node v) {
            current_blossom[v] = _b;
        });
    }
    for (auto [b, e] : blossom->subblossoms) {
        if (b->label == even) {
            calc_best_edges(b);
        }
    }
}

void GabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        Blossom* _b = b;
        b->for_nodes([this, _b] (NetworKit::node v) {
            current_blossom[v] = _b;
        });
    }
}

void GabowMaximumMatching::adjust_by_delta(MaximumMatching::edgeweight delta) {
    graph.forNodes([this, delta] (NetworKit::node v) {
        Blossom* v_blossom = this->get_blossom(v);
        if (v_blossom->label == even) {
            U[v] -= delta;
        } else if (v_blossom->label == odd) {
            U[v] += delta;
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

MaximumMatching::edgeweight GabowMaximumMatching::calc_delta1() {
    // min u_i : i - even vertex

    MaximumMatching::edgeweight res = std::numeric_limits<MaximumMatching::edgeweight>::max();
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (this->get_blossom(v)->label == even) {
            res = std::min(res, this->U[v]);
        }
    });
    return res;
}

MaximumMatching::edgeweight GabowMaximumMatching::calc_delta2() {
    // min pi_ij : i - even vertex, j - free vertex

    MaximumMatching::edgeweight res = std::numeric_limits<MaximumMatching::edgeweight>::max();
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (get_blossom(v)->label != free) return;
        auto slack = edge_slack(best_edge[v].id);
        res = std::min(res, slack);
    });
    return res;
}

MaximumMatching::edgeweight GabowMaximumMatching::calc_delta3() {
    // min pi_ij / 2 : i,j - even vertices in different blossoms

    MaximumMatching::edgeweight res = std::numeric_limits<MaximumMatching::edgeweight>::max();
    for (auto b : blossoms) {
        if (b->label != even) continue;
        auto slack = edge_slack(get_data(b)->best_edge.id);
        res = std::min(res, slack / 2);
    }
    return res;
}

MaximumMatching::edgeweight GabowMaximumMatching::calc_delta4() {
    // min z_k / 2 : B_k - odd blossom 

    MaximumMatching::edgeweight res = std::numeric_limits<MaximumMatching::edgeweight>::max();
    for (Blossom* b : blossoms) {
        if (b->label == odd && !b->is_trivial()) {
            res = std::min(res, b->z / 2);
        }
    }

    return res;
}

void GabowMaximumMatching::find_delta2_useful_edges() {
    graph.forNodes([this] (NetworKit::node v) {
        if (get_blossom(v)->label != free) return;
        if (edge_slack(best_edge[v].id) == 0) 
            edge_queue.push(best_edge[v]);
    });
}

void GabowMaximumMatching::find_delta3_useful_edges() {
    for (auto b : blossoms) {
        if (b->label != even) continue;
        auto edge = get_data(b)->best_edge;
        if (edge_slack(edge.id) == 0) 
            edge_queue.push(edge);
    }
}

std::vector<GabowMaximumMatching::Blossom*> GabowMaximumMatching::get_odd_blossoms_to_expand() {
    std::vector<Blossom*> to_expand;
    for (auto b : blossoms) {
        if (b->label == odd && b->z == 0.0 && !b->is_trivial()) {
            to_expand.push_back(b);
        }
    }
    return to_expand;
}

GabowMaximumMatching::Blossom* GabowMaximumMatching::get_blossom(NetworKit::node vertex) {
    return current_blossom[vertex];
}

void GabowMaximumMatching::calc_best_edges(Blossom* b) {
    auto data = get_data(b);

    data->best_edges.clear();
    data->best_edge = no_edge;
    auto best_slack = edge_slack(data->best_edge.id);

    b->for_nodes([this, b, data, &best_slack] (NetworKit::node u) {
        graph.forEdgesOf(u, [this, b, data, &best_slack] 
                (NetworKit::node u, NetworKit::node v, NetworKit::node id) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) return;

            auto slack = edge_slack(id);
            if (v_blossom->label == even) {
                if (data->best_edges.find(v_blossom) == data->best_edges.end() 
                    || slack < edge_slack(data->best_edges[v_blossom].id)) {
                    data->best_edges[v_blossom] = {u, v, id};
                    if (slack < best_slack) {
                        best_slack = slack;
                        data->best_edge = {u, v, id};
                    }
                }
            } else {
                if (slack < edge_slack(best_edge[v].id)) {
                    best_edge[v] = {u, v, id};
                }
            }
        });
    });
}

MaximumMatching::edgeweight GabowMaximumMatching::edge_slack(NetworKit::edgeid edge) {
    if (edge == NetworKit::none) return std::numeric_limits<MaximumMatching::edgeweight>::max();
    auto [u, v, w] = graph_edges[edge];
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);
    return U[u] + U[v] - w + (u_blossom == v_blossom ? u_blossom->z : 0);
}

void GabowMaximumMatching::check_consistency() {
    std::cerr << "Current vertices: \n";
    graph.forNodes([this] (NetworKit::node v) {
        if (this->matched_vertex[v] != NetworKit::none)
            std::cerr << v << " " << this->matched_vertex[v] << "  ";
        else std::cerr << v << " none  ";
        this->get_blossom(v)->short_print(); std::cerr << std::endl;
    });

    std::cerr << "Current weights: \n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << "U[" << v << "] = " << U[v] << std::endl;
    });

    for (auto b : blossoms) if (!b->is_trivial()) {
        b->short_print(); std::cerr << ": " << b->z << std::endl; 
    }

    std::cerr << "Current best edges: \n";
    graph.forNodes([this] (NetworKit::node v) {
        if (get_blossom(v)->label != even) {
            std::cerr << "best_edge[" << v << "] = (" 
                << best_edge[v].u << ", " << best_edge[v].v << ") : " 
                << edge_slack(best_edge[v].id) << std::endl;
        }
    });
    for (auto b : blossoms) {
        if (b->label != even) continue;
        std::cerr << "Best edges of "; b->short_print(); std::cerr << std::endl;
        for (auto [c, e] : get_data(b)->best_edges) {
            c->short_print(); 
            std::cerr << " : (" << e.u << ", " << e.v << ") : " << edge_slack(e.id) << std::endl;
        }
        auto edge = get_data(b)->best_edge;
        std::cerr << "Best one: (" << edge.u << ", " << edge.v << ") : " << edge_slack(edge.id) << std::endl;
    }
}

} /* namespace Koala */
