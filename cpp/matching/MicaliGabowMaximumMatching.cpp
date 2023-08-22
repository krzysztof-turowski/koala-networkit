#include <matching/MaximumMatching.hpp>

namespace Koala {

MicaliGabowMaximumMatching::MicaliGabowMaximumMatching(NetworKit::Graph &graph) : 
        BlossomMaximumMatching(graph),
        nodes_refs(graph.upperNodeIdBound()),
        Ueven(graph.upperNodeIdBound()),
        Uodd(graph.upperNodeIdBound()),
        Zeven(graph.upperNodeIdBound()),
        Zodd(graph.upperNodeIdBound()),
        good_edges(graph.upperEdgeIdBound()),
        even_edges(graph.upperEdgeIdBound() + graph.upperNodeIdBound()) { 

    NetworKit::edgeweight max_weight = std::numeric_limits<NetworKit::edgeweight>::min();
    graph.forEdges(
        [&max_weight] (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) { 
            max_weight = std::max(weight, max_weight);
        }
    );
    Ufree = std::vector<NetworKit::edgeweight>(graph.upperNodeIdBound(), max_weight / 2.);

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
    Uodd.clear();
    Ueven.clear();
    good_edges.clear();
    Zeven.clear();
    Zodd.clear();

    for (auto blossom : blossoms) {
        blossom->label = is_exposed(blossom) ? even : free;
        blossom->backtrack_edge = no_edge;
        
        if (blossom->label == free) {
            auto even_edge_group = even_edges.new_group(true);
            blossom->for_nodes([this, even_edge_group] (NetworKit::node v) {
                even_edges.append_dummy(dummy_edge_id(v), even_edge_group);
            });
            get_data(blossom)->even_edges = even_edge_group;
        }
        if (blossom->label == even) {
            if (!blossom->is_trivial()) {
                Zeven.insert(blossom->initial_base, blossom->z);
            }
            blossom->for_nodes([this] (NetworKit::node u) {
                Ueven.insert(u, Ufree[u]);
            });
        }
    }

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

    Zeven.for_elements([this, &to_expand] (NetworKit::node base, NetworKit::edgeweight dual_weight) {
        Blossom* b = get_blossom(base);
        b->z = dual_weight;
        if (dual_weight == 0)
            to_expand.push_back(b);
    });
    Zodd.for_elements([this] (NetworKit::node base, NetworKit::edgeweight dual_weight) {
        Blossom* b = get_blossom(base);
        b->z = dual_weight;
    });
    Ueven.for_elements([this] (NetworKit::node v, NetworKit::edgeweight dual) {
        Ufree[v] = dual;
    });
    Uodd.for_elements([this] (NetworKit::node v, NetworKit::edgeweight dual) {
        Ufree[v] = dual;
    });
    
    for (auto b : to_expand) expand_even_blossom(b);
}

void MicaliGabowMaximumMatching::initialize_substage() {}

bool MicaliGabowMaximumMatching::has_useful_edges() {
    clear_not_good_edges();
    return (!good_edges.empty() && good_edges.find_min().second == 0) || 
           (even_edges.has_active_elements() && even_edges.find_min().second == 0);
}

MicaliGabowMaximumMatching::EdgeInfo MicaliGabowMaximumMatching::get_useful_edge() {
    if (!good_edges.empty() && good_edges.find_min().second == 0) {
        auto id = good_edges.find_min().first;
        auto [u, v, w] = graph_edges[id];
        good_edges.remove_min();

        return {u, v, id};
    }

    if (even_edges.has_active_elements() && even_edges.find_min().second == 0) {
        auto id = even_edges.find_min().first;
        auto [u, v, w] = graph_edges[id];
        even_edges.remove(id);

        return {u, v, id};
    }

    return no_edge;
}

void MicaliGabowMaximumMatching::label_odd(Blossom* b) {
    b->for_nodes([this] (NetworKit::node v) {
        Uodd.insert(v, Ufree[v]);
    });
    if (!b->is_trivial()) { 
        Zodd.insert(b->initial_base, b->z);
    }
    even_edges.change_status(get_data(b)->even_edges, false);
}

void MicaliGabowMaximumMatching::label_even(Blossom* b) {
    b->for_nodes([this, b] (NetworKit::node u) {
        Ueven.insert(u, Ufree[u]);
    });
    if (!b->is_trivial()) {
        Zeven.insert(b->initial_base, b->z);
    }
    even_edges.delete_group(get_data(b)->even_edges);
    get_data(b)->even_edges = nullptr;
    scan_edges(b);
}

void MicaliGabowMaximumMatching::scan_edges(Blossom* b) {
    b->for_nodes([this, b] (NetworKit::node u) {
        graph.forEdgesOf(u, [this, b] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            auto v_blossom = get_blossom(v);
            if (v_blossom == b) return;
            auto slack = edge_slack(id);
            if (v_blossom->label == even) {
                good_edges.insert(id, slack/2);
            } else {
                even_edges.insert_before(
                    id, slack, dummy_edge_id(v), get_data(v_blossom)->even_edges);
            }
        });
    });
}

void MicaliGabowMaximumMatching::clear_not_good_edges() {
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
    ConcatenableQueue<Blossom*, NetworKit::node, NetworKit::node> nodes(new_blossom);
    for (auto [b, e] : new_blossom->sub_blossoms) {
        nodes.concat(std::move(get_data(b)->nodes), new_blossom);

        if (b->label == even && !b->is_trivial()) {
            b->z = Zeven.current_priority(b->initial_base);
            Zeven.remove(b->initial_base);
        } else if (b->label == odd) {
            if (!b->is_trivial()) {
                b->z = Zodd.current_priority(b->initial_base);
                Zodd.remove(b->initial_base);
            }
            b->for_nodes([this] (NetworKit::node v) {
                auto u = Uodd.current_priority(v);
                Uodd.remove(v);
                Ueven.insert(v, u);
            });
        }
    }
    for (auto [b, e] : new_blossom->sub_blossoms) {
        if (b->label == odd) {
            even_edges.delete_group(get_data(b)->even_edges);
            get_data(b)->even_edges = nullptr;
            scan_edges(b);
        }
    }
    new_blossom->data = new MicaliGabowBlossomData(std::move(nodes), nullptr);
    Zeven.insert(new_blossom->initial_base, 0);
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

    Zodd.remove(blossom->initial_base);

    for (auto [b, e] : blossom->sub_blossoms) {
        auto [nodes_b, nodes_rest] = remaining_nodes->split(nodes_refs[b->last_node], b, blossom);
        remaining_nodes = nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b); delete nodes_b;
        
        auto [edges_b, edges_rest] = even_edges.split_group(remaining_edges, dummy_edge_id(b->last_node));
        remaining_edges = edges_rest;
        if (b->label == even) {
            even_edges.delete_group(edges_b);
        } else {
            get_data(b)->even_edges = edges_b;
            even_edges.change_status(edges_b, b->label == free);
        }

        if (b->label == even) {
            b->for_nodes([this] (NetworKit::node v) {
                auto dual = Uodd.current_priority(v);
                Uodd.remove(v);
                Ueven.insert(v, dual);
            });
        } else if (b->label == free) {
            b->for_nodes([this] (NetworKit::node v) {
                auto dual = Uodd.current_priority(v);
                Uodd.remove(v);
                Ufree[v] = dual;
            });
        }

        if (b->label == even && !b->is_trivial()) {
            Zeven.insert(b->initial_base, b->z);
        }
        if (b->label == odd && !b->is_trivial()) {
            Zodd.insert(b->initial_base, b->z);
        }
    }

    for (auto [b, e] : blossom->sub_blossoms) {
        if (b->label == even) {
            scan_edges(b);
        }
    }
}

void MicaliGabowMaximumMatching::handle_even_blossom_expansion(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    auto remaining_nodes = &get_data(blossom)->nodes;
    
    for (auto [b, e] : blossom->sub_blossoms) {
        auto [nodes_b, nodes_rest] = remaining_nodes->split(nodes_refs[b->last_node], b, blossom);
        remaining_nodes = nodes_rest;
        get_data(b)->nodes = std::move(*nodes_b); delete nodes_b;
    }
}

void MicaliGabowMaximumMatching::adjust_by_delta(NetworKit::edgeweight delta) {
    Ueven.decrease_all_priorities(delta);
    Uodd.decrease_all_priorities(-delta);
    good_edges.decrease_all_priorities(delta);
    Zeven.decrease_all_priorities(-2*delta);
    Zodd.decrease_all_priorities(2*delta);
    even_edges.decrease_all_priorities(delta);
}

NetworKit::edgeweight MicaliGabowMaximumMatching::calc_delta1() {
    // min u_i : i - even vertex
    return Ueven.empty() ? std::numeric_limits<NetworKit::edgeweight>::max()
        : Ueven.find_min().second;
}

NetworKit::edgeweight MicaliGabowMaximumMatching::calc_delta2() {
    // min pi_ij : i - even vertex, j - free vertex

    return even_edges.has_active_elements() ? even_edges.find_min().second
        : std::numeric_limits<NetworKit::edgeweight>::max();
}

NetworKit::edgeweight MicaliGabowMaximumMatching::calc_delta3() {
    // min pi_ij / 2 : i,j - even vertices in different blossoms
    
    clear_not_good_edges();
    return good_edges.empty() ? std::numeric_limits<NetworKit::edgeweight>::max()
        : good_edges.find_min().second;
}

NetworKit::edgeweight MicaliGabowMaximumMatching::calc_delta4() {
    // min z_k / 2 : B_k - odd blossom 
    return Zodd.empty() ? std::numeric_limits<NetworKit::edgeweight>::max()
        : (Zodd.find_min().second / 2);
}

void MicaliGabowMaximumMatching::find_delta2_useful_edges() {}

void MicaliGabowMaximumMatching::find_delta3_useful_edges() {}

std::vector<MicaliGabowMaximumMatching::Blossom*> 
MicaliGabowMaximumMatching::get_odd_blossoms_to_expand() {    
    std::vector<Blossom*> to_expand;
    Zodd.for_elements_until([this, &to_expand] (NetworKit::node base, NetworKit::edgeweight dual) {
        if (dual == 0) {
            to_expand.push_back(get_blossom(base));
            return false;
        }
        return true;
    });

    return to_expand;
}

MicaliGabowMaximumMatching::Blossom* 
MicaliGabowMaximumMatching::get_blossom(NetworKit::node vertex) {
    return nodes_refs[vertex]->find_queue()->head;
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
        std::cerr << ": " << blossom_dual(blossom) << std::endl;
    }
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << ": " << U(v) << std::endl;
    });
    std::cerr << "Zeven:\n";
    Zeven.for_elements([this] (NetworKit::node base, NetworKit::edgeweight dual_weight) {
        std::cerr << base << " : " << dual_weight << std::endl;
    });
    std::cerr << "Zodd:\n";
    Zodd.for_elements([this] (NetworKit::node base, NetworKit::edgeweight dual_weight) {
        std::cerr << base << " : " << dual_weight << std::endl;
    });
    std::cerr << "Good edges:\n";
    good_edges.for_elements([this] (NetworKit::edgeid id, NetworKit::edgeweight var) {
        auto [u, v, w] = graph_edges[id];
        std::cerr << "(" << u << ", " << v << ") : " << var << std::endl;
    });
    std::cerr << "Even edges:\n";
    for (auto blossom : blossoms) {
        if (blossom->label == even) continue;        
        blossom->nodes_print(); std::cerr << std::endl;
        EvenEdgeGroup group = get_data(blossom)->even_edges;
        std::cerr << (group->active ? "active" : "non-active") << std::endl;
        if (group->empty()) {
            std::cerr << "Empty\n";
        } else {
            auto [id, pi] = group->find_min();
            auto [u, v, w] = graph_edges[id];
                std::cerr << "min : (" << u << ", " << v << ") : " << pi << std::endl;        
        }
        even_edges.for_each_in_group(group, 
            [this] (NetworKit::edgeid id, NetworKit::edgeweight pi) {
            if (id < graph.upperEdgeIdBound()) {
                auto [u, v, w] = graph_edges[id];
                std::cerr << "(" << u << ", " << v << ") : " << pi << std::endl;        
            } else {
                // std::cerr << "dummy " << id - graph.upperEdgeIdBound() << std::endl;
            }
        });
    }
    std::cerr << "Active minima:\n";
    even_edges.for_group_minima([this] (NetworKit::edgeid id, NetworKit::edgeweight pi) {
        if (id < graph.upperEdgeIdBound()) {
            auto [u, v, w] = graph_edges[id];
            std::cerr << "(" << u << ", " << v << ") : " << pi << std::endl;        
        } else {
            // std::cerr << "dummy " << id - graph.upperEdgeIdBound() << std::endl;
        }
    });
}

NetworKit::edgeweight MicaliGabowMaximumMatching::U(NetworKit::node v) {
    auto b = get_blossom(v);
    switch (b->label) {
        case free: return Ufree[v];
        case even: return Ueven.current_priority(v);
        case odd:  return Uodd.current_priority(v);
    }
    return 0;
}

NetworKit::edgeweight MicaliGabowMaximumMatching::edge_slack(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);
    return U(u) + U(v) - w + (u_blossom == v_blossom ? blossom_dual(u_blossom) : 0);
}

NetworKit::edgeweight MicaliGabowMaximumMatching::blossom_dual(Blossom* b) {
    switch (b->label) {
        case free: return b->z;
        case even: return Zeven.current_priority(b->initial_base);
        case odd:  return Zodd.current_priority(b->initial_base);
    }
    return 0;
}

NetworKit::edgeid MicaliGabowMaximumMatching::dummy_edge_id(NetworKit::node node) {
    return graph.upperEdgeIdBound() + node;
}

} /* namespace Koala */
