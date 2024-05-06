#include <matching/MaximumMatching.hpp>

namespace Koala {

MicaliGabowMaximumMatching::MicaliGabowMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        nodes_refs(graph.upperNodeIdBound()),
        y_even(graph.upperNodeIdBound()),
        y_odd(graph.upperNodeIdBound()),
        y_free(graph.upperNodeIdBound(), max_weight / 2),
        z_even(graph.upperNodeIdBound()),
        z_odd(graph.upperNodeIdBound()),
        good_edges(graph.upperEdgeIdBound()),
        even_edges(graph.upperEdgeIdBound() + graph.upperNodeIdBound()) { 

    for (auto b : trivial_blossom) {
        ConcatenableQueue<Blossom*, NetworKit::node, NetworKit::node> nodes(b);
        nodes_refs[b->base] = nodes.append(b->base, 0);
        b->data = new MicaliGabowBlossomData(std::move(nodes), nullptr);
    }
}

MicaliGabowMaximumMatching::MicaliGabowBlossomData* 
MicaliGabowMaximumMatching::get_data(Blossom* b) {
    return static_cast<MicaliGabowBlossomData*>(b->data);
}

void MicaliGabowMaximumMatching::initialize_stage() {
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
            blossom->for_nodes([this, even_edge_group] (NetworKit::node v) {
                // Insert dummy edges for each node in blossom order
                // This is needed when an odd blossom is expanded to allow splitting
                even_edges.append_dummy(dummy_edge_id(v), even_edge_group);
            });
            get_data(blossom)->even_edges = even_edge_group;
        }

        if (blossom->label == even) {
            if (!blossom->is_trivial()) {
                // Track dual weight for non trivial even blossoms
                z_even.insert(blossom->initial_base, blossom->z);
            }
            // Track weights for even vertices
            blossom->for_nodes([this] (NetworKit::node u) {
                y_even.insert(u, y_free[u]);
            });
        }
    }

    // Scan outgoing edges for all initial even blossoms
    for (auto blossom : blossoms) {
        if (blossom->label == even) {
            scan_edges(blossom);
        } 
    }
}

void MicaliGabowMaximumMatching::finish_stage() {
    // Expand even blossoms with B_z = 0
    std::vector<Blossom*> to_expand;

    for (auto blossom : blossoms) {
        if (blossom->label != even) {
            even_edges.delete_group(get_data(blossom)->even_edges);
            get_data(blossom)->even_edges = nullptr;
        }
    }

    // Retrieve final weights for blossoms and vertices at the end of the stage
    z_even.for_elements([this, &to_expand] (NetworKit::node base, MaximumWeightMatching::edgeweight dual_weight) {
        Blossom* b = get_blossom(base);
        b->z = dual_weight;
        if (dual_weight == 0)
            to_expand.push_back(b);
    });
    z_odd.for_elements([this] (NetworKit::node base, MaximumWeightMatching::edgeweight dual_weight) {
        Blossom* b = get_blossom(base);
        b->z = dual_weight;
    });
    y_even.for_elements([this] (NetworKit::node v, MaximumWeightMatching::edgeweight dual) {
        y_free[v] = dual;
    });
    y_odd.for_elements([this] (NetworKit::node v, MaximumWeightMatching::edgeweight dual) {
        y_free[v] = dual;
    });
    
    for (auto b : to_expand) expand_even_blossom(b);
}

void MicaliGabowMaximumMatching::initialize_substage() {}

bool MicaliGabowMaximumMatching::has_useful_edges() {
    clear_not_good_edges();
    
    return (!good_edges.empty() && good_edges.find_min().second == 0) || 
           (even_edges.has_active_elements() && even_edges.find_min().second == 0);
}

MicaliGabowMaximumMatching::Edge MicaliGabowMaximumMatching::get_useful_edge() {
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
        auto id = even_edges.find_min().first;
        auto [u, v, w] = graph_edges[id];
        even_edges.remove(id);

        return {u, v, id};
    }

    return no_edge;
}

void MicaliGabowMaximumMatching::handle_grow(Blossom* odd_blossom, Blossom* even_blossom) {
    // Begin tracking dual weights for newly odd vertices
    odd_blossom->for_nodes([this] (NetworKit::node v) {
        y_odd.insert(v, y_free[v]);
    });
    if (!odd_blossom->is_trivial()) { 
        // Track dual weight for the newly odd blossom
        z_odd.insert(odd_blossom->initial_base, odd_blossom->z);
    }
    
    // Even edges from the odd blossom are no longer affected by dual adjustments
    // Deactivate the corresponding group
    even_edges.change_status(get_data(odd_blossom)->even_edges, false);

    // Begin tracking dual weights for newly even vertices
    even_blossom->for_nodes([this, even_blossom] (NetworKit::node u) {
        y_even.insert(u, y_free[u]);
    });
    if (!even_blossom->is_trivial()) {
        // Track dual weight for the newly even blossom
        z_even.insert(even_blossom->initial_base, even_blossom->z);
    }
    // Delete the group as the blossom is now even
    even_edges.delete_group(get_data(even_blossom)->even_edges);
    get_data(even_blossom)->even_edges = nullptr;

    // Scan edges from the even blossom
    scan_edges(even_blossom);
}

void MicaliGabowMaximumMatching::scan_edges(Blossom* b) {
    b->for_nodes([this, b] (NetworKit::node u) {
        graph.forEdgesOf(u, [this, b] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) return;

            auto edge_slack = slack(id);

            if (v_blossom->label == even) {
                // A good edge is found. Begin tracking it
                good_edges.insert(id, edge_slack);
            } else {
                // An even edge is found. Add it to the group corresponding to the non-even blossom
                // Preserve the order by inserting before a dummy node
                even_edges.insert_before(
                    id, edge_slack, dummy_edge_id(v), get_data(v_blossom)->even_edges);
            }
        });
    });
}

void MicaliGabowMaximumMatching::clear_not_good_edges() {
    // Lazily clear edges that are no longer good
    while (!good_edges.empty()) {
        auto [edge, _] = good_edges.find_min();
        if (!is_good(edge)) 
            good_edges.remove_min();
        else break;
    }
}

bool MicaliGabowMaximumMatching::is_good(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    return get_blossom(u) != get_blossom(v);
}

void MicaliGabowMaximumMatching::handle_new_blossom(Blossom* new_blossom) {
    // Create a list of nodes for the new blossom to track it's vertices
    ConcatenableQueue<Blossom*, NetworKit::node, NetworKit::node> nodes(new_blossom);

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

            b->for_nodes([this] (NetworKit::node v) {
                // Change which priority queue handles the dual weight for a previously odd vertex
                auto u = y_odd.current_priority(v);
                y_odd.remove(v);
                y_even.insert(v, u);
            });
        }
    }

    for (auto [b, e] : new_blossom->subblossoms) {
        if (b->label == odd) {
            // Delete the even edge group as the blossom is no longer odd
            even_edges.delete_group(get_data(b)->even_edges);
            get_data(b)->even_edges = nullptr;

            // Scan the vertices that become even for the first time to find good and even edges
            scan_edges(b);
        }
    }

    new_blossom->data = new MicaliGabowBlossomData(std::move(nodes), nullptr);

    // Track the dual weight for the new blossom
    z_even.insert(new_blossom->initial_base, 0);
}

void MicaliGabowMaximumMatching::handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) {    
    auto data = get_data(blossom);

    auto nodes = std::move(data->nodes);
    auto [nodesA, nodesB] = nodes.split(nodes_refs[subblossom->last_node], blossom, blossom);
    nodesB->concat(std::move(*nodesA), blossom);
    data->nodes = std::move(*nodesB);
    delete nodesA;
    delete nodesB;

    if (data->even_edges != nullptr) {
        even_edges.shift_group(data->even_edges, dummy_edge_id(subblossom->last_node));
    }
}

void MicaliGabowMaximumMatching::handle_odd_blossom_expansion(Blossom* blossom) {
    auto remaining_nodes = &get_data(blossom)->nodes;
    auto remaining_edges = get_data(blossom)->even_edges;

    // Notice we don't remove blossom from z_odd as it's been done in get_odd_blossoms_to_expand()

    for (auto [b, e] : blossom->subblossoms) {
        // Split the node list of the expanded blossoms
        // Assign each subblossom a list 
        auto [nodes_b, nodes_rest] = remaining_nodes->split(nodes_refs[b->last_node], b, blossom);
        remaining_nodes = nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b); 
        delete nodes_b;
        
        // Split even edge groups of the expanded blossom
        // Use the dummy edges to achieve the task
        auto [edges_b, edges_rest] = even_edges.split_group(remaining_edges, dummy_edge_id(b->last_node));
        remaining_edges = edges_rest;
        if (b->label == even) {
            // Delete the group if it now corresponds to an even blossom
            even_edges.delete_group(edges_b);
        } else {
            // Assign the group 
            get_data(b)->even_edges = edges_b;

            // Wether a group is now active is de
            even_edges.change_status(edges_b, b->label == free);
        }

        // Retrieve the current weight for vertices that stop being odd and change how they're tracked
        if (b->label == even) {
            b->for_nodes([this] (NetworKit::node v) {
                auto dual = y_odd.current_priority(v);
                y_odd.remove(v);
                y_even.insert(v, dual);
            });
        } else if (b->label == free) {
            b->for_nodes([this] (NetworKit::node v) {
                auto dual = y_odd.current_priority(v);
                y_odd.remove(v);
                y_free[v] = dual;
            });
        }

        // Start tracking dual weight for non-free subblossoms
        if (b->label == even && !b->is_trivial()) {
            z_even.insert(b->initial_base, b->z);
        }
        if (b->label == odd && !b->is_trivial()) {
            z_odd.insert(b->initial_base, b->z);
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

void MicaliGabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    auto remaining_nodes = &get_data(blossom)->nodes;
    
    for (auto [b, e] : blossom->subblossoms) {
        auto [nodes_b, nodes_rest] = remaining_nodes->split(nodes_refs[b->last_node], b, blossom);
        remaining_nodes = nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b); delete nodes_b;
    }
}

void MicaliGabowMaximumMatching::adjust_by_delta(MaximumWeightMatching::edgeweight delta) {
    // Update dual weight for even and odd vertices
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

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::calc_delta1() {
    // Find the even vertex with minimum dual weight
    // min u_i : i - even vertex
    return y_even.empty() ? std::numeric_limits<MaximumWeightMatching::edgeweight >::max()
        : y_even.find_min().second;
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::calc_delta2() {
    // Find the even edge with the smallest slack
    // min pi_ij : i - even vertex, j - free vertex

    return even_edges.has_active_elements() ? even_edges.find_min().second
        : std::numeric_limits<MaximumWeightMatching::edgeweight >::max();
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::calc_delta3() {
    // Find the good edge with smallest slack
    // min pi_ij / 2 : i,j - even vertices in different blossoms
    
    // Some edges might no longer be good and have to be removed
    clear_not_good_edges();
    return good_edges.empty() ? std::numeric_limits<MaximumWeightMatching::edgeweight>::max()
        : good_edges.find_min().second / 2;
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::calc_delta4() {
    // Find the odd blossom with minimum dual weight
    // min z_k / 2 : B_k - odd blossom 

    return z_odd.empty() ? std::numeric_limits<MaximumWeightMatching::edgeweight>::max()
        : z_odd.find_min().second / 2;
}

void MicaliGabowMaximumMatching::find_delta2_useful_edges() {}

void MicaliGabowMaximumMatching::find_delta3_useful_edges() {}

std::vector<MicaliGabowMaximumMatching::Blossom*> 
MicaliGabowMaximumMatching::get_odd_blossoms_to_expand() {    
    std::vector<Blossom*> to_expand;
    while (!z_odd.empty() && z_odd.find_min().second == 0) {
        to_expand.push_back(get_blossom(z_odd.find_min().first));
        z_odd.remove_min();
    }

    return to_expand;
}

MicaliGabowMaximumMatching::Blossom* 
MicaliGabowMaximumMatching::get_blossom(NetworKit::node vertex) {
    return nodes_refs[vertex]->find_queue()->head;
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::y(NetworKit::node v) {
    // Find the current dual weight for a vertex depending on it's label
    auto b = get_blossom(v);
    switch (b->label) {
        case free: return y_free[v];
        case even: return y_even.current_priority(v);
        case odd:  return y_odd.current_priority(v);
    }
    return 0;
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::slack(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);
    return y(u) + y(v) - w + (u_blossom == v_blossom ? z(u_blossom) : 0);
}

MaximumWeightMatching::edgeweight  MicaliGabowMaximumMatching::z(Blossom* b) {
    // Find the dual weight for a blossom depending on it's label
    switch (b->label) {
        case free: return b->z;
        case even: return z_even.current_priority(b->initial_base);
        case odd:  return z_odd.current_priority(b->initial_base);
    }
    return 0;
}

NetworKit::edgeid MicaliGabowMaximumMatching::dummy_edge_id(NetworKit::node node) {
    return graph.upperEdgeIdBound() + node;
}

void MicaliGabowMaximumMatching::check_consistency() {
    std::cerr << "Current blossoms:\n";
    for (auto blossom : blossoms) {
        get_data(blossom)->nodes.check_consistency();
        if (blossom->is_trivial()) continue;
        blossom->short_print(); std::cerr << std::endl;
        get_data(blossom)->nodes.for_each([this, blossom] (NetworKit::node v, NetworKit::node _) {
            std::cerr << v << " ";
            if (get_blossom(v) != blossom) {
                std::cerr << "get_blossom(" << v << ") has wrong value\n";
                exit(1);
            }
        });
        std::cerr << std::endl;
    }
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << ": ";
        get_blossom(v)->nodes_print();
        std::cerr << std::endl;
    });
    std::cerr << "Current weights:\n";
    for (auto blossom : blossoms) {
        if (blossom->is_trivial()) continue;
        blossom->nodes_print();
        std::cerr << ": " << z(blossom) << std::endl;
    }
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << ": " << y(v) << std::endl;
    });
    std::cerr << "z_even:\n";
    z_even.for_elements([this] (NetworKit::node base, MaximumWeightMatching::edgeweight dual_weight) {
        std::cerr << base << " : " << dual_weight << std::endl;
    });
    std::cerr << "z_odd:\n";
    z_odd.for_elements([this] (NetworKit::node base, MaximumWeightMatching::edgeweight dual_weight) {
        std::cerr << base << " : " << dual_weight << std::endl;
    });
    std::cerr << "y_even:\n";
    y_even.for_elements([this] (NetworKit::node v, MaximumWeightMatching::edgeweight dual_weight) {
        std::cerr << v << " : " << dual_weight << std::endl;
    });
    std::cerr << "y_odd:\n";
    y_odd.for_elements([this] (NetworKit::node v, MaximumWeightMatching::edgeweight dual_weight) {
        std::cerr << v << " : " << dual_weight << std::endl;
    });
    std::cerr << "Good edges: " << good_edges.heap.size() << "\n";
    good_edges.for_elements([this] (NetworKit::edgeid id, MaximumWeightMatching::edgeweight var) {
        auto [u, v, w] = graph_edges[id];
        std::cerr << "(" << u << ", " << v << ") : " << var << std::endl;
    });
    std::cerr << "Even edges:\n";
    for (auto blossom : blossoms) {
        if (blossom->label == even) continue;        
        
        EvenEdgeGroup group = get_data(blossom)->even_edges;
        std::cerr << "Group for ";
        blossom->nodes_print(); std::cerr << " - ";
        std::cerr << (group->active ? "active" : "non-active") << std::endl;
        if (group->empty()) {
            std::cerr << "Empty\n";
        } else {
            auto [id, pi] = group->find_min();
            auto [u, v, w] = graph_edges[id];
                std::cerr << "min : (" << u << ", " << v << ") : " << pi << std::endl;        
        }
        even_edges.for_each_in_group(group, 
            [this] (NetworKit::edgeid id, MaximumWeightMatching::edgeweight pi) {
            if (id < graph.upperEdgeIdBound()) {
                auto [u, v, w] = graph_edges[id];
                std::cerr << "(" << u << ", " << v << ") : " << pi << std::endl;        
            } else {
                // std::cerr << "dummy " << id - graph.upperEdgeIdBound() << std::endl;
            }
        });
    }
    std::cerr << "Active minima:\n";
    even_edges.for_group_minima([this] (NetworKit::edgeid id, MaximumWeightMatching::edgeweight pi) {
        if (id < graph.upperEdgeIdBound()) {
            auto [u, v, w] = graph_edges[id];
            std::cerr << "(" << u << ", " << v << ") : " << pi << std::endl;        
        } else {
            // std::cerr << "dummy " << id - graph.upperEdgeIdBound() << std::endl;
        }
    });
}

} /* namespace Koala */
