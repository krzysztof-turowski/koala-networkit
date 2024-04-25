#include <matching/MaximumMatching.hpp>

namespace Koala {

GabowMaximumMatching::GabowMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        current_blossom(graph.upperNodeIdBound(), nullptr),
        U(graph.upperNodeIdBound(), 0),
        best_edge(graph.upperNodeIdBound(), no_edge),
        sorted_neighbours(graph.upperNodeIdBound()) { 
    
    graph.sortEdges();
    MaximumMatching::edgeweight max_weight = std::numeric_limits<MaximumMatching::edgeweight>::min();
    graph.forEdges([this, &max_weight] 
            (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w, NetworKit::edgeid id) { 
        max_weight = std::max(static_cast<MaximumMatching::edgeweight>(w), max_weight);
        sorted_neighbours[u].emplace_back(v, id);
        sorted_neighbours[v].emplace_back(u, id);
    });
    U = std::vector<MaximumMatching::edgeweight>(graph.upperNodeIdBound(), max_weight);
    graph.forNodes([this] (NetworKit::node vertex) {
        current_blossom[vertex] = trivial_blossom[vertex];
        std::sort(sorted_neighbours[vertex].begin(), sorted_neighbours[vertex].end());
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

void GabowMaximumMatching::handle_grow(Blossom* odd_blossom, Blossom* even_blossom) {
    calc_best_edges(even_blossom);
}

void GabowMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    for (auto [b, edge] : new_blossom->subblossoms) { 
        b->for_nodes([this, new_blossom] (NetworKit::node v) {
            current_blossom[v] = new_blossom;
        });
    }
    GabowBlossomData* data = new GabowBlossomData();
    new_blossom->data = data;
    data->best_edge = no_edge;
    auto best_slack = edge_slack(NetworKit::none);

    #if DEBUG_LOGGING
    std::cerr << "Creating a list for ";
    new_blossom->nodes_print();
    std::cerr << std::endl;
    #endif

    for (auto [b, e] : new_blossom->subblossoms) {
        if (b->label == odd) 
            calc_best_edges(b);

        #if DEBUG_LOGGING
        std::cerr << "Current list:\n";
        for (auto [u, v, id] : data->best_edges) {
            std::cerr << "(" << u << ", " << v << ") : " << edge_slack(id) << std::endl;
        }
        std::cerr << "Merging list of ";
        b->nodes_print();
        std::cerr << std::endl;
        #endif

        auto b_data = get_data(b);
        auto edge_it = data->best_edges.begin();

        for (auto edge : b_data->best_edges) {
            auto [u, v, id] = edge;
            if (get_blossom(v) == new_blossom) continue;

            auto slack = edge_slack(id);

            #if DEBUG_LOGGING
            std::cerr << "check (" << u << ", " << v << ") : " << slack << std::endl;
            #endif
            
            while (edge_it != data->best_edges.end() && v > edge_it->v) edge_it ++;

            if (edge_it == data->best_edges.end() || v < edge_it->v) {
                #if DEBUG_LOGGING
                std::cerr << "insert" << std::endl;
                #endif

                data->best_edges.insert(edge_it, {u, v, id});
            } else if (slack < edge_slack(edge_it->id)) {
                #if DEBUG_LOGGING
                std::cerr << "replace (" << edge_it->u << ", " << edge_it->v << ") : " << edge_slack(edge_it->id) << std::endl;
                #endif

                *edge_it = {u, v, id};
            }

            if (slack < best_slack) {
                best_slack = slack;
                data->best_edge = edge;
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

    #if DEBUG_LOGGING
    // std::cerr << "Calculating best edges for ";
    // b->short_print();
    // std::cerr << std::endl;
    #endif

    b->for_nodes([this, b, data, &best_slack] (NetworKit::node u) {
        auto edge_it = data->best_edges.begin();

        // Iterate in sorted order for efficiency
        for (auto [v, id] : sorted_neighbours[u]) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) continue;

            auto slack = edge_slack(id);
            if (v_blossom->label == even) {
                // Find first edge (u', v') where v' <= v
                while (edge_it != data->best_edges.end() && v > edge_it->v) edge_it ++;

                if (edge_it == data->best_edges.end() || v < edge_it->v) {
                    // No edge (u', v)
                    data->best_edges.insert(edge_it, {u, v, id});
                } else if (slack < edge_slack(edge_it->id)) {
                    // Only insert if smaller slack
                    *edge_it = {u, v, id};
                }

                // Update best slack
                if (slack < best_slack) {
                    best_slack = slack;
                    data->best_edge = {u, v, id};
                }
            } else {
                // Update best edge for an outer vertex
                if (slack < edge_slack(best_edge[v].id)) {
                    best_edge[v] = {u, v, id};
                }
            }
        }
    });

    #if DEBUG_LOGGING
    // for (auto [u, v, id] : data->best_edges) {
    //     std::cerr << "(" << u << ", " << v << ") : " << edge_slack(id) << std::endl;
    // }
    // auto [u, v, id] = data->best_edge;
    // std::cerr << "Best : " << "(" << u << ", " << v << ") : " << edge_slack(id) << std::endl;
    #endif
}

MaximumMatching::edgeweight GabowMaximumMatching::edge_slack(NetworKit::edgeid edge) {
    if (edge == NetworKit::none) return std::numeric_limits<MaximumMatching::edgeweight>::max();
    auto [u, v, w] = graph_edges[edge];
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);

    #if DEBUG_LOGGING
    // std::cerr << "edge_slack(" << u << ", " << v << ") = " << U[u] << " + " << U[v] << " - " << w << " + ?" << std::endl;
    #endif

    return U[u] + U[v] - w + (u_blossom == v_blossom ? u_blossom->z : 0);
}

void GabowMaximumMatching::check_consistency() {
    std::cerr << "Current vertices: \n";
    graph.forNodes([this] (NetworKit::node v) {
        if (this->matched_vertex[v] != NetworKit::none)
            std::cerr << v << " " << this->matched_vertex[v] << " : ";
        else std::cerr << v << " - : ";
        this->get_blossom(v)->nodes_print(); 
        std::cerr << std::endl;
    });

    std::cerr << "Current weights: \n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << ": " << U[v] << std::endl;
    });

    for (auto b : blossoms) if (!b->is_trivial()) {
        b->short_print(); 
        std::cerr << " : " << b->z << std::endl; 
    }

    std::cerr << "Current best edges: \n";
    graph.forNodes([this] (NetworKit::node v) {
        if (get_blossom(v)->label != even) {
            std::cerr << "best_edge[" << v << "] = " 
                << edge_to_string(best_edge[v]) << " : " 
                << edge_slack(best_edge[v].id) << std::endl;
        }
    });
    for (auto b : blossoms) {
        if (b->label != even) continue;
        std::cerr << "Best edges of "; b->short_print(); std::cerr << std::endl;
        for (auto e : get_data(b)->best_edges) {
            std::cerr << edge_to_string(e) << " : " << edge_slack(e.id) << std::endl;
        }
        auto edge = get_data(b)->best_edge;
        std::cerr << "Best one: " << edge_to_string(edge) << " : " << edge_slack(edge.id) << std::endl;
    }
}

} /* namespace Koala */
