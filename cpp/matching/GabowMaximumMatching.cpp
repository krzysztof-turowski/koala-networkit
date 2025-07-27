#include <matching/MaximumMatching.hpp>

namespace Koala {

GabowMaximumMatching::GabowMaximumMatching(
    NetworKit::Graph &graph, bool perfect, InitializationStrategy initialization) :
        BlossomMaximumMatching(graph, perfect, initialization),
        current_blossom(graph.upperNodeIdBound(), nullptr),
        best_edge(graph.upperNodeIdBound(), no_edge) {
    // Sort the edges for efficient best edge list contstruction
    graph.sortEdges();

    graph.forNodes([this] (NetworKit::node vertex) {
        current_blossom[vertex] = trivial_blossom[vertex];
        trivial_blossom[vertex]->data = new GabowBlossomData();
    });
}

GabowMaximumMatching::GabowBlossomData* GabowMaximumMatching::get_data(Blossom* b) {
    return static_cast<GabowBlossomData*>(b->data);
}

void GabowMaximumMatching::initialize_stage() {
    edge_queue = {};

    graph.forNodes([this] (NetworKit::node v) {
        best_edge[v] = no_edge;
    });

    for (auto blossom : blossoms) {
        // Start search in all exposed blossoms
        blossom->label = is_exposed(blossom) ? even : free;

        // Clear the list of best edges
        auto data = get_data(blossom);
        data->best_edges.clear();
        blossom->backtrack_edge = data->best_edge = no_edge;
    }
    for (auto blossom : blossoms) {
        if (blossom->label == even) {
            scan_edges(blossom);
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

GabowMaximumMatching::Edge GabowMaximumMatching::get_useful_edge() {
    auto edge = edge_queue.front();
    edge_queue.pop();
    return edge;
}

void GabowMaximumMatching::handle_grow(Blossom* odd_blossom, Blossom* even_blossom) {
    // Scan the edges for newly even vertices
    scan_edges(even_blossom);
}

void GabowMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    GabowBlossomData* data = new GabowBlossomData();
    new_blossom->data = data;
    data->best_edge = no_edge;
    auto best_slack = slack(NetworKit::none);

    for (auto [b, e] : new_blossom->subblossoms)
        for (auto v : b->nodes)
            current_blossom[v] = new_blossom;

    for (auto [b, e] : new_blossom->subblossoms) {
        if (b->label == odd)
            scan_edges(b);

        auto b_data = get_data(b);
        auto edge_it = data->best_edges.begin();

        // Merge subblossom lists of best edges into the list for new edge
        for (auto edge : b_data->best_edges) {
            auto [u, v, id] = edge;
            if (get_blossom(v) == new_blossom) continue;

            auto edge_slack = slack(id);

            // Find first edge (u', v') where v' >= v
            while (edge_it != data->best_edges.end() && v > edge_it->v) edge_it++;

            if (edge_it == data->best_edges.end() || v < edge_it->v) {
                // No edge (u', v)
                data->best_edges.insert(edge_it, {u, v, id});
            } else if (edge_slack < slack(edge_it->id)) {
                // Replace existing (u', v) if smaller slack
                *edge_it = {u, v, id};
            }

            // Update best slack
            if (edge_slack < best_slack) {
                best_slack = edge_slack;
                data->best_edge = edge;
            }
        }

        b_data->best_edges.clear();
        b_data->best_edge = no_edge;
    }
}

void GabowMaximumMatching::handle_nodes_split(Blossom* blossom) {
    for (auto [b, e] : blossom->subblossoms) {
        for (auto v : b->nodes) {
            current_blossom[v] = b;
        }
    }
}

void GabowMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    // Scan edges for newly even vertices
    for (auto [b, e] : blossom->subblossoms) {
        if (b->label == even) {
            scan_edges(b);
        }
    }
}

void GabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {}

void GabowMaximumMatching::adjust_by_delta(MaximumWeightMatching::weight delta) {
    graph.forNodes([this, delta] (NetworKit::node v) {
        Blossom* v_blossom = this->get_blossom(v);
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

MaximumWeightMatching::weight GabowMaximumMatching::calc_delta1() {
    // min u_i : i - even vertex

    MaximumWeightMatching::weight res = infinite_weight;
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (get_blossom(v)->label == even) {
            res = std::min(res, y[v]);
        }
    });
    return res;
}

MaximumWeightMatching::weight GabowMaximumMatching::calc_delta2() {
    // min pi_ij : i - even vertex, j - free vertex

    MaximumWeightMatching::weight res = infinite_weight;
    graph.forNodes([this, &res] (NetworKit::node v) {
        if (get_blossom(v)->label != free) return;
        auto edge_slack = slack(best_edge[v].id);
        res = std::min(res, edge_slack);
    });
    return res;
}

MaximumWeightMatching::weight GabowMaximumMatching::calc_delta3() {
    // min pi_ij / 2 : i,j - even vertices in different blossoms

    MaximumWeightMatching::weight res = infinite_weight;
    for (auto b : blossoms) {
        if (b->label != even) continue;
        auto edge_slack = slack(get_data(b)->best_edge.id);
        res = std::min(res, edge_slack / 2);
    }
    return res;
}

MaximumWeightMatching::weight GabowMaximumMatching::calc_delta4() {
    // min z_k / 2 : B_k - odd blossom

    MaximumWeightMatching::weight res = infinite_weight;
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
        if (slack(best_edge[v].id) == 0)
            edge_queue.push(best_edge[v]);
    });
}

void GabowMaximumMatching::find_delta3_useful_edges() {
    for (auto b : blossoms) {
        if (b->label != even) continue;
        auto edge = get_data(b)->best_edge;
        if (slack(edge.id) == 0)
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

void GabowMaximumMatching::scan_edges(Blossom* b) {
    auto data = get_data(b);

    data->best_edges.clear();
    data->best_edge = no_edge;
    auto best_slack = slack(data->best_edge.id);

    for (auto u : b->nodes) {
        auto edge_it = data->best_edges.begin();

        // Iterate in sorted order for efficiency
        graph.forEdgesOf(u, [this, b, data, &best_slack, &edge_it]
                (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) return;

            auto edge_slack = slack(id);
            if (v_blossom->label == even) {
                // Find first edge (u', v') where v' >= v
                while (edge_it != data->best_edges.end() && v > edge_it->v) edge_it++;

                if (edge_it == data->best_edges.end() || v < edge_it->v) {
                    // No edge (u', v)
                    data->best_edges.insert(edge_it, {u, v, id});
                } else if (edge_slack < slack(edge_it->id)) {
                    // Replace existing (u', v) if smaller slack
                    *edge_it = {u, v, id};
                }

                // Update best slack
                if (edge_slack < best_slack) {
                    best_slack = edge_slack;
                    data->best_edge = {u, v, id};
                }
            } else {
                // Update best edge for an outer vertex
                if (edge_slack < slack(best_edge[v].id)) {
                    best_edge[v] = {u, v, id};
                }
            }
        });
    }
}

MaximumWeightMatching::weight GabowMaximumMatching::slack(NetworKit::edgeid edge) {
    if (edge == NetworKit::none) return infinite_weight;
    auto [u, v, w] = graph_edges[edge];
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);

    return y[u] + y[v] - w + (u_blossom == v_blossom ? u_blossom->z : 0);
}

} /* namespace Koala */
