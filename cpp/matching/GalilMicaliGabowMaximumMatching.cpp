#include <matching/MaximumMatching.hpp>

namespace Koala {

GalilMicaliGabowMaximumMatching::GalilMicaliGabowMaximumMatching(
    NetworKit::Graph &graph, bool perfect, InitializationStrategy initialization):
        BlossomMaximumMatching(graph, perfect, initialization),
        nodes_refs(graph.upperNodeIdBound()),
        y_even(graph.upperNodeIdBound()),
        y_odd(graph.upperNodeIdBound()),
        y_free(y),
        z_even(graph.upperNodeIdBound()),
        z_odd(graph.upperNodeIdBound()),
        good_edges(graph.upperEdgeIdBound()),
        even_edges(graph.upperNodeIdBound(), NetworKit::none, NetworKit::none, infinite_weight) {
    for (auto b : trivial_blossom) {
        BlossomNodeList nodes(b);
        nodes_refs[b->base] = nodes.append(b->base, 0);
        b->data = new GalilMicaliGabowBlossomData(std::move(nodes), nullptr);
    }
}

GalilMicaliGabowMaximumMatching::GalilMicaliGabowBlossomData*
GalilMicaliGabowMaximumMatching::get_data(Blossom* b) {
    return static_cast<GalilMicaliGabowBlossomData*>(b->data);
}

void GalilMicaliGabowMaximumMatching::initialize_stage() {
    edge_queue = {};
    y_odd.clear();
    y_even.clear();
    good_edges.clear();
    z_even.clear();
    z_odd.clear();

    for (auto blossom : blossoms) {
        // Mark all exposed blossoms as even
        blossom->label = is_exposed(blossom) ? even : free;
        blossom->backtrack_edge = no_edge;

        if (blossom->label == free) {
            // Create edge groups for all non even blossoms
            auto even_edge_group = even_edges.new_group(true);
            for (auto v : blossom->nodes) {
                // Insert dummy edges for each node in blossom order
                // This is needed when an odd blossom is expanded to allow splitting
                even_edges.append(v, NetworKit::none, infinite_weight, even_edge_group);
            }
            get_data(blossom)->even_edge_group = even_edge_group;
        }

        if (blossom->label == even) {
            if (!blossom->is_trivial()) {
                // Track dual weight for non trivial even blossoms
                z_even.insert(blossom->initial_base, blossom, blossom->z);
            }
            // Track weights for even vertices
            for (auto v : blossom->nodes) {
                y_even.insert(v, v, y_free[v]);
            }
        }
    }

    // Scan outgoing edges for all initial even blossoms
    for (auto blossom : blossoms) {
        if (blossom->label == even) {
            scan_edges(blossom);
        }
    }
}

void GalilMicaliGabowMaximumMatching::finish_stage() {
    // Expand even blossoms with B_z = 0
    std::vector<Blossom*> to_expand;

    for (auto blossom : blossoms) {
        if (blossom->label != even) {
            even_edges.delete_group(get_data(blossom)->even_edge_group);
            get_data(blossom)->even_edge_group = nullptr;
        }

        // Retrieve final weights for blossoms and vertices at the end of the stage
        if (blossom->label == even) {
            blossom->z = z_even.current_priority(blossom->initial_base);
            for (auto v : blossom->nodes) {
                y_free[v] = y_even.current_priority(v);
            }
            // Check if blossom has to be expanded
            if (blossom->z == 0)
                to_expand.push_back(blossom);
        } else if (blossom->label == odd) {
            blossom->z = z_odd.current_priority(blossom->initial_base);
            for (auto v : blossom->nodes) {
                y_free[v] = y_odd.current_priority(v);
            }
        }
    }

    for (auto b : to_expand) expand_even_blossom(b);
}

void GalilMicaliGabowMaximumMatching::initialize_substage() {}

bool GalilMicaliGabowMaximumMatching::has_useful_edges() {
    clear_not_good_edges();

    return (!good_edges.empty() && good_edges.find_min().second == 0) ||
           (even_edges.has_active_elements() && even_edges.find_min().second == 0);
}

GalilMicaliGabowMaximumMatching::Edge GalilMicaliGabowMaximumMatching::get_useful_edge() {
    // Find tight edges

    // Check if there are good edges and the one with minimum slack is tight
    if (!good_edges.empty() && good_edges.find_min().second == 0) {
        auto id = good_edges.find_min().first;
        auto [u, v, w] = graph_edges[id];
        good_edges.remove_min();

        return {u, v, id};
    }

    // If there are no tight good edges try to do the same for even edges
    if (even_edges.has_active_elements() && even_edges.find_min().second == 0) {
        auto [vertex, id, slack] = even_edges.find_min_element();
        auto [u, v, w] = graph_edges[id];
        // No need to remove the minimum as the grow step will change the blossom to odd and
        // deactivate the group

        return {u, v, id};
    }

    return no_edge;
}

void GalilMicaliGabowMaximumMatching::handle_grow(Blossom* odd_blossom, Blossom* even_blossom) {
    // Begin tracking dual weights for newly odd vertices
    for (auto v : odd_blossom->nodes) {
        y_odd.insert(v, v, y_free[v]);
    }
    if (!odd_blossom->is_trivial()) {
        // Track dual weight for the newly odd blossom
        z_odd.insert(odd_blossom->initial_base, odd_blossom, odd_blossom->z);
    }

    // Even edges from the odd blossom are no longer affected by dual adjustments
    // Deactivate the corresponding group
    even_edges.change_status(get_data(odd_blossom)->even_edge_group, false);

    // Begin tracking dual weights for newly even vertices
    for (auto v : even_blossom->nodes) {
        y_even.insert(v, v, y_free[v]);
    }
    if (!even_blossom->is_trivial()) {
        // Track dual weight for the newly even blossom
        z_even.insert(even_blossom->initial_base, even_blossom, even_blossom->z);
    }
    // Delete the group as the blossom is now even
    even_edges.delete_group(get_data(even_blossom)->even_edge_group);
    get_data(even_blossom)->even_edge_group = nullptr;

    // Scan edges from the even blossom
    scan_edges(even_blossom);
}

void GalilMicaliGabowMaximumMatching::scan_edges(Blossom* b) {
    for (auto u : b->nodes) {
        graph.forEdgesOf(u, [this, b] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) return;

            auto edge_slack = slack(id);

            if (v_blossom->label == even) {
                // A good edge is found. Begin tracking it
                good_edges.insert(id, id, edge_slack);
            } else {
                // An even edge is found. Update priority for v in the even edge group
                even_edges.decrease_priority(v, id, edge_slack,
                    get_data(v_blossom)->even_edge_group);
            }
        });
    }
}

void GalilMicaliGabowMaximumMatching::clear_not_good_edges() {
    // Lazily clear edges that are no longer good
    while (!good_edges.empty()) {
        auto [edge, _] = good_edges.find_min();
        if (!is_good(edge))
            good_edges.remove_min();
        else
            break;
    }
}

bool GalilMicaliGabowMaximumMatching::is_good(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    return get_blossom(u) != get_blossom(v);
}

void GalilMicaliGabowMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    // Create a list of nodes for the new blossom to track it's vertices
    BlossomNodeList nodes(new_blossom);

    for (auto [b, e] : new_blossom->subblossoms) {
        // Concatenate the lists of all the subblossoms
        nodes.concat(std::move(get_data(b)->nodes), new_blossom);

        if (b->label == even && !b->is_trivial()) {
            // Retrieve the current weight for the subblossom
            b->z = z_even.current_priority(b->initial_base);
            z_even.remove(b->initial_base);
        } else if (b->label == odd) {
            if (!b->is_trivial()) {
                // Retrieve the current weight for the subblossom
                b->z = z_odd.current_priority(b->initial_base);
                z_odd.remove(b->initial_base);
            }

            for (auto v : b->nodes) {
                // Change which priority queue handles the dual weight for a previously odd vertex
                auto u = y_odd.current_priority(v);
                y_odd.remove(v);
                y_even.insert(v, v, u);
            }
        }
    }

    for (auto [b, e] : new_blossom->subblossoms) {
        if (b->label == odd) {
            // Delete the even edge group as the blossom is no longer odd
            even_edges.delete_group(get_data(b)->even_edge_group);
            get_data(b)->even_edge_group = nullptr;

            // Scan the vertices that become even for the first time to find good and even edges
            scan_edges(b);
        }
    }

    new_blossom->data = new GalilMicaliGabowBlossomData(std::move(nodes), nullptr);

    // Track the dual weight for the new blossom
    z_even.insert(new_blossom->initial_base, new_blossom, 0);
}

void GalilMicaliGabowMaximumMatching::handle_subblossom_shift(
        Blossom* blossom, Blossom* subblossom) {
    auto data = get_data(blossom);

    // Change the order of nodes in the blossom node
    BlossomNodeList nodes(std::move(data->nodes));
    auto [nodesA, nodesB] = nodes.split(nodes_refs[subblossom->last_node], blossom, blossom);
    nodesB->concat(std::move(*nodesA), blossom);
    data->nodes = std::move(*nodesB);
    delete nodesA;
    delete nodesB;

    if (data->even_edge_group != nullptr) {
        // Change the order of nodes in blossom's even edge group
        even_edges.shift_group(data->even_edge_group, subblossom->last_node);
    }
}

void GalilMicaliGabowMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    auto data = get_data(blossom);
    auto remaining_edges = data->even_edge_group;
    // Notice we don't remove blossom from z_odd as it's been done in get_odd_blossoms_to_expand()

    for (auto [b, e] : blossom->subblossoms) {
        // Split the node list of the expanded blossoms
        // Assign each subblossom a list
        auto [nodes_b, nodes_rest] = data->nodes.split(nodes_refs[b->last_node], b, blossom);
        get_data(blossom)->nodes = std::move(*nodes_rest);
        delete nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b);
        delete nodes_b;

        // Split even edge groups of the expanded blossom
        // Use the dummy edges to achieve the task
        auto [edges_b, edges_rest] = even_edges.split_group(remaining_edges, b->last_node);
        remaining_edges = edges_rest;
        if (b->label == even) {
            // Delete the group if it now corresponds to an even blossom
            even_edges.delete_group(edges_b);
        } else {
            // Assign the group
            get_data(b)->even_edge_group = edges_b;

            // Wether a group is now active is de
            even_edges.change_status(edges_b, b->label == free);
        }

        // Retrieve the current weight for vertices that stop being odd and change which queue
        // they're tracked by
        if (b->label == even) {
            for (auto v : b->nodes) {
                auto dual = y_odd.current_priority(v);
                y_odd.remove(v);
                y_even.insert(v, v, dual);
            }
        } else if (b->label == free) {
            for (auto v : b->nodes) {
                auto dual = y_odd.current_priority(v);
                y_odd.remove(v);
                y_free[v] = dual;
            }
        }

        // Start tracking dual weight for non-free subblossoms
        if (b->label == even && !b->is_trivial()) {
            z_even.insert(b->initial_base, b, b->z);
        }
        if (b->label == odd && !b->is_trivial()) {
            z_odd.insert(b->initial_base, b, b->z);
        }
    }

    even_edges.delete_group(remaining_edges);
    get_data(blossom)->even_edge_group = nullptr;

    // Scan the vertices that become even for the first time to find good and even edges
    // Do this only after all the other data is consistent
    for (auto [b, e] : blossom->subblossoms) {
        if (b->label == even) {
            scan_edges(b);
        }
    }
}

void GalilMicaliGabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    auto data = get_data(blossom);

    for (auto [b, e] : blossom->subblossoms) {
        auto [nodes_b, nodes_rest] = data->nodes.split(nodes_refs[b->last_node], b, blossom);
        data->nodes = std::move(*nodes_rest);
        delete nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b);
        delete nodes_b;
    }
}

void GalilMicaliGabowMaximumMatching::adjust_by_delta(MaximumWeightMatching::weight delta) {
    // Update dual weights for even and odd vertices
    y_even.decrease_all_priorities(delta);
    y_odd.decrease_all_priorities(-delta);

    // Update dual weights for even and odd blossoms
    z_even.decrease_all_priorities(-2 * delta);
    z_odd.decrease_all_priorities(2 * delta);

    // Update slack for good edges
    good_edges.decrease_all_priorities(2 * delta);

    // Update slack for even edges
    even_edges.decrease_all_priorities(delta);
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::calc_delta1() {
    // Find the even vertex with minimum dual weight
    // min u_i : i - even vertex
    return y_even.empty() ? infinite_weight : y_even.find_min().second;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::calc_delta2() {
    // Find the even edge with the smallest slack
    // min pi_ij : i - even vertex, j - free vertex
    return even_edges.has_active_elements() ? even_edges.find_min().second : infinite_weight;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::calc_delta3() {
    // Find the good edge with smallest slack
    // min pi_ij / 2 : i,j - even vertices in different blossoms
    // Some edges might no longer be good and have to be removed
    clear_not_good_edges();
    return good_edges.empty() ? infinite_weight : good_edges.find_min().second / 2;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::calc_delta4() {
    // Find the odd blossom with minimum dual weight
    // min z_k / 2 : B_k - odd blossom
    return z_odd.empty() ? infinite_weight : z_odd.find_min().second / 2;
}

void GalilMicaliGabowMaximumMatching::find_delta2_useful_edges() {}

void GalilMicaliGabowMaximumMatching::find_delta3_useful_edges() {}

std::vector<GalilMicaliGabowMaximumMatching::Blossom*>
GalilMicaliGabowMaximumMatching::get_odd_blossoms_to_expand() {
    std::vector<Blossom*> to_expand;
    while (!z_odd.empty() && z_odd.find_min().second == 0) {
        to_expand.push_back(z_odd.find_min().first);
        z_odd.remove_min();
    }

    return to_expand;
}

GalilMicaliGabowMaximumMatching::Blossom*
GalilMicaliGabowMaximumMatching::get_blossom(NetworKit::node vertex) {
    return nodes_refs[vertex]->find_queue()->root_id;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::current_y(NetworKit::node v) {
    // Find the current dual weight for a vertex depending on it's label
    auto b = get_blossom(v);
    switch (b->label) {
        case free: return y_free[v];
        case even: return y_even.current_priority(v);
        case odd:  return y_odd.current_priority(v);
    }
    return 0;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::slack(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    return current_y(u) + current_y(v) - w;
}

} /* namespace Koala */
