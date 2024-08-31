#include <matching/MaximumMatching.hpp>

namespace Koala {

GalilMicaliGabowMaximumMatching::GalilMicaliGabowMaximumMatching(
    NetworKit::Graph &graph, bool perfect, InitializationStrategy initialization):
        BlossomMaximumMatching(graph, perfect, initialization),
        y_even(graph.upperNodeIdBound()),
        z_even(graph.upperNodeIdBound()),
        z_odd(graph.upperNodeIdBound()),
        good_edges(graph.upperEdgeIdBound()),
        even_edges(graph.upperNodeIdBound(), NetworKit::none, NetworKit::none, infinite_weight) {
    for (auto b : trivial_blossom) {
        auto group = even_edges.new_group(b);
        even_edges.append(b->base, NetworKit::none, infinite_weight, group);
        b->data = new GalilMicaliGabowBlossomData(group);
    }
}

GalilMicaliGabowMaximumMatching::GalilMicaliGabowBlossomData*
GalilMicaliGabowMaximumMatching::get_data(Blossom* b) {
    return static_cast<GalilMicaliGabowBlossomData*>(b->data);
}

void GalilMicaliGabowMaximumMatching::initialize_stage() {
    edge_queue = {};
    y_even.clear();
    good_edges.clear();
    z_even.clear();
    z_odd.clear();
    even_edges.reset();

    for (auto blossom : blossoms) {
        // Mark all exposed blossoms as even
        blossom->label = is_exposed(blossom) ? even : free;
        blossom->backtrack_edge = no_edge;
        auto group = get_data(blossom)->node_group;
        even_edges.reset_group(group);

        if (blossom->label == free) {
            // Create edge groups for all non even blossoms
            even_edges.change_status(group, true);
        }

        if (blossom->label == even) {
            if (!blossom->is_trivial()) {
                // Track dual weight for non trivial even blossoms
                z_even.insert(blossom->initial_base, blossom, blossom->z);
            }
            // Track weights for even vertices
            for (auto v : blossom->nodes) {
                y_even.insert(v, v, y[v]);
            }
        } else {
            auto group = get_data(blossom)->node_group;
            for (auto v : blossom->nodes) {
                even_edges.set_element_priority(v, y[v], group);
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
        // Retrieve final weights for blossoms and vertices at the end of the stage
        if (blossom->label == even) {
            blossom->z = z_even.current_priority(blossom->initial_base);
            for (auto v : blossom->nodes) {
                y[v] = y_even.current_priority(v);
            }
            // Check if blossom has to be expanded
            if (blossom->z == 0)
                to_expand.push_back(blossom);
        } else if (blossom->label == odd) {
            blossom->z = z_odd.current_priority(blossom->initial_base);
        }

        if (blossom->label != even) {
            auto group = get_data(blossom)->node_group;
            for (auto v : blossom->nodes) {
                y[v] = even_edges.current_element_priority(v, group);
            }
        }
    }

    for (auto b : to_expand) expand_even_blossom(b);
}

void GalilMicaliGabowMaximumMatching::initialize_substage() {
    // std::cerr << "y:\n";
    // graph.forNodes([this] (NetworKit::node v) {
    //     std::cerr << v << ": " << current_y(v) << std::endl;
    // });
}

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
    // for (auto v : odd_blossom->nodes) {
    //     even_edges.set_element_priority(v, y[v], get_data(odd_blossom)->node_group);
    // }
    if (!odd_blossom->is_trivial()) {
        // Track dual weight for the newly odd blossom
        z_odd.insert(odd_blossom->initial_base, odd_blossom, odd_blossom->z);
    }

    // Even edges from the odd blossom are no longer affected by dual adjustments
    // Deactivate the corresponding group
    even_edges.change_status(get_data(odd_blossom)->node_group, false);

    // Begin tracking dual weights for newly even vertices
    for (auto v : even_blossom->nodes) {
        auto y_v = even_edges.current_element_priority(v, get_data(even_blossom)->node_group);
        y_even.insert(v, v, y_v);
    }
    if (!even_blossom->is_trivial()) {
        // Track dual weight for the newly even blossom
        z_even.insert(even_blossom->initial_base, even_blossom, even_blossom->z);
    }
    // Kill the group as the blossom is now even
    even_edges.change_status(get_data(even_blossom)->node_group, false);

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
                    get_data(v_blossom)->node_group);
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
    // Create a group for the new blossom to track it's vertices
    auto group = even_edges.new_group(new_blossom);

    for (auto [b, e] : new_blossom->subblossoms) {
        auto group_b = get_data(b)->node_group;

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
                auto u = even_edges.current_element_priority(v, group_b);
                y_even.insert(v, v, u);
            }
        }

        if (b->label == odd) {
            even_edges.change_status(group_b, false);
        }

        // Concatenate the lists of all the subblossoms
        group = even_edges.concat_groups(group, group_b, new_blossom);
        get_data(b)->node_group = nullptr;
    }

    for (auto [b, e] : new_blossom->subblossoms) {
        if (b->label == odd) {
            // Scan the vertices that become even for the first time to find good and even edges
            scan_edges(b);
        }
    }

    new_blossom->data = new GalilMicaliGabowBlossomData(group);

    // Track the dual weight for the new blossom
    z_even.insert(new_blossom->initial_base, new_blossom, 0);
}

void GalilMicaliGabowMaximumMatching::handle_nodes_split(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    auto group = get_data(blossom)->node_group;
    even_edges.change_status(group, false);

    for (auto [b, e] : blossom->subblossoms) {
        // Split the node list of the expanded blossoms
        // Assign each subblossom a list
        auto [group_b, rest] = even_edges.split_group(group, b->last_node, b, blossom);
        group = rest;
        get_data(b)->node_group = group_b;
    }

    delete group;
    get_data(blossom)->node_group = nullptr;
}

void GalilMicaliGabowMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    // Notice we don't remove blossom from z_odd as it's been done in get_odd_blossoms_to_expand()
    // The node lists have been updated in handle_nodes_split

    for (auto [b, e] : blossom->subblossoms) {
        auto group_b = get_data(b)->node_group;

        // Split even edge groups of the expanded blossom
        // Use the dummy edges to achieve the task
        if (b->label == even) {
            // Kill the group if it now corresponds to an even blossom
            even_edges.change_status(group_b, false);
        } else {
            // Whether a group is now active is dependent on the label
            even_edges.change_status(group_b, b->label == free);
        }

        // Retrieve the current weight for vertices that stop being odd and change which queue
        // they're tracked by
        if (b->label == even) {
            for (auto v : b->nodes) {
                auto dual = even_edges.current_element_priority(v, group_b);
                y_even.insert(v, v, dual);
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

    // Scan the vertices that become even for the first time to find good and even edges
    // Do this only after all the other data is consistent
    for (auto [b, e] : blossom->subblossoms) {
        if (b->label == even) {
            scan_edges(b);
        }
    }
}

void GalilMicaliGabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {}

void GalilMicaliGabowMaximumMatching::adjust_by_delta(MaximumWeightMatching::weight delta) {
    // Update dual weights for even and odd vertices
    y_even.decrease_all_priorities(delta);

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
    return even_edges.find_group(vertex);
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::current_y(NetworKit::node v) {
    // Find the current dual weight for a vertex depending on it's label
    auto b = get_blossom(v);
    switch (b->label) {
        case even:
            return y_even.current_priority(v);
        case free:
        case odd:
            return even_edges.current_element_priority(v, get_data(b)->node_group);
    }
    return 0;
}

MaximumWeightMatching::weight GalilMicaliGabowMaximumMatching::slack(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    return current_y(u) + current_y(v) - w;
}

} /* namespace Koala */
