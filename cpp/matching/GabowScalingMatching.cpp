#include <matching/MaximumMatching.hpp>

namespace Koala {

std::string node_to_str(NetworKit::node v) {
    return v == NetworKit::none ? "-" : std::to_string(v);
}

NetworKit::Graph reduce_to_MWPM(const NetworKit::Graph& graph) {
    // Reduce the graph to and instance of Maximum Weight Perfect Matching.
    // The reduced graph consists of two copies of the original graph were corresponding 
    // vertices are connected with a zero weight edges.
    auto original_graph_size = graph.upperNodeIdBound();
    NetworKit::Graph G(2 * original_graph_size, true, false, true);
    graph.forEdges(
        [&G, original_graph_size] (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
            G.addEdge(u, v, w);
            G.addEdge(u + original_graph_size, v + original_graph_size, w);
        }
    );
    graph.forNodes([&G, original_graph_size] (NetworKit::node v) {
        G.addEdge(v, v + original_graph_size, 0);
    });
    return G;
}

GabowScalingMatching::GabowScalingMatching(NetworKit::Graph &graph):
    MaximumWeightMatching(graph), 
    reducedGraph(reduce_to_MWPM(graph)), 
    current_y(reducedGraph.upperNodeIdBound()),
    trivial_blossom(reducedGraph.upperNodeIdBound(), nullptr),
    union_find(reducedGraph.upperNodeIdBound()),
    matched_vertex(reducedGraph.upperNodeIdBound()),
    matched_edge(reducedGraph.upperNodeIdBound()),
    split_find_min(reducedGraph.upperNodeIdBound(), 1000000000, no_edge, 
        reducedGraph.numberOfNodes(), reducedGraph.numberOfEdges()),
    y0(reducedGraph.upperNodeIdBound()),
    Delta(reducedGraph.upperNodeIdBound()),
    t_shell(reducedGraph.upperNodeIdBound()),
    current_blossom(reducedGraph.upperNodeIdBound(), nullptr),
    current_shell(reducedGraph.upperNodeIdBound()),
    search_shell(reducedGraph.upperNodeIdBound()),
    graph_edges(reducedGraph.upperEdgeIdBound()),
    actual_to_contracted(reducedGraph.upperEdgeIdBound()),
    edge_in_matching(reducedGraph.upperEdgeIdBound()),
    vertex_path(reducedGraph.upperNodeIdBound(), nullptr) { 
        reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
            graph_edges[e] = {u, v};
        });
    }

void GabowScalingMatching::run() {
    std::vector<MaximumWeightMatching::intweight> w(reducedGraph.upperEdgeIdBound());
    reducedGraph.forEdges(
        [&w] (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid e) {
            // Multiply all weights by 2 to ensure they're even
            w[e] = 2 * static_cast<MaximumWeightMatching::intweight>(ew);
        }
    );

    auto [M, _, __] = scale(w);

    graph.forNodes([this, &M] (NetworKit::node v) {
        // Only include matched vertices if they're from the original graph and weren't added 
        // in the reduction
        matching[v] = M[v] < graph.upperNodeIdBound() ? M[v] : NetworKit::none;
    });

    hasRun = true;
}

bool allWeightsZero(const std::vector<int>& w) {
    for (auto wi : w) if (wi != 0) return false;
    return true;
}

std::vector<NetworKit::node> defaultMatching(const NetworKit::Graph& original_graph) {
    auto original_graph_size = original_graph.upperNodeIdBound();
    std::vector<NetworKit::node> M(2 * original_graph_size, NetworKit::none);
    original_graph.forNodes([&M, original_graph_size] (NetworKit::node v) {
        M[v] = v + original_graph_size;
        M[v + original_graph_size] = v;
    });
    return M;
}

std::tuple<std::vector<NetworKit::node>, std::vector<int>, GabowScalingMatching::OldBlossom*>
GabowScalingMatching::scale(const std::vector<int>& w) {
    if (allWeightsZero(w)) {
        // Base case - all weights are zero
        // Return only an outer blossom and a zero weight for all vertices
        // The perfect matching consists of edges joining corresponding vertices in the reduced graph
        OldBlossom* T = new OldBlossom;
        T->size = reducedGraph.numberOfNodes();
        T->z = 0;
        T->heavy_path_parent = nullptr;
        T->heavy_child = nullptr;
        T->dissolved = false;
        reducedGraph.forNodes([T] (NetworKit::node v) {
            T->nodes.push_back(v);
        });
        return std::make_tuple(
            defaultMatching(graph),
            std::vector<int>(reducedGraph.upperNodeIdBound(), 0),
            T
        );
    }

    // Scale down weights
    std::vector<int> w1(w.size());
    for (int i = 0; i < w.size(); ++ i) 
        w1[i] = 2 * (w[i] / 4); // Keep all weights even

    // Recursive call to scale for reduced weights
    auto [_, y1, T] = scale(w1);

    current_w = w;
    old_root = T;
    T->for_blossoms([this] (OldBlossom* B) { 
        // Use weights from recursive call
        B->z = 2 * B->z; 
    });
    reducedGraph.forNodes([this, &y1] (NetworKit::node v) {
        // Use weights from recursive call
        current_y[v] = 2 * y1[v] + 1;
        matched_vertex[v] = NetworKit::none;
        matched_edge[v] = NetworKit::none;
    });
    // Reset matching
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        edge_in_matching[e] = false;
    });
    edges_in_matching = 0;
    
    // Find matching by dissolving old blossoms
    match(T);

    // Convert found blossoms into the representation for old blossoms
    auto new_T = turn_current_blossoms_into_old(T->shell_blossoms);

    return std::make_tuple(matched_vertex, current_y, new_T);
}

void GabowScalingMatching::match(OldBlossom* T) {
    create_trivial_blossoms(T);
    
    heavy_path_decomposition(T, 0);
}

void GabowScalingMatching::delete_blossom(Blossom* B, OldBlossom* S) {
    // Delete a blossom from the shell
    S->shell_blossoms.erase(B->shell_blossoms_it);
}

void GabowScalingMatching::add_blossom(Blossom* B, OldBlossom* S) {
    // Add a blossom to the shell
    S->shell_blossoms.push_back(B);
    B->shell_blossoms_it = std::prev(S->shell_blossoms.end());
}

void GabowScalingMatching::expand_blossom(Blossom* B, OldBlossom* S) {
    if (B->is_trivial()) return;
    
    // Add subblossom of the expanded blossom to their parent's shell
    for (auto [b, e] : B->subblossoms) {
        b->parent = nullptr;
        add_blossom(b, S);
    }

    // Delete the expanded blossom from the shell
    delete_blossom(B, S);
    delete B;
}

void GabowScalingMatching::swap_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v] = graph_edges[edge];
    
    if (edge_in_matching[edge]) {
        edge_in_matching[edge] = false;
        edges_in_matching --;
        if (matched_vertex[u] == v)
            matched_vertex[u] = matched_edge[u] = NetworKit::none;
        if (matched_vertex[v] == u)
            matched_vertex[v] = matched_edge[v] = NetworKit::none;
    } else {
        edge_in_matching[edge] = true;
        edges_in_matching ++;
        matched_vertex[u] = v;
        matched_vertex[v] = u;
        matched_edge[u] = matched_edge[v] = edge;
    }
}

void GabowScalingMatching::set_edge_in_matching(NetworKit::edgeid edge) {
    if (edge_in_matching[edge]) {
        check_edge_in_matching(edge);
        return;
    }

    auto [u, v] = graph_edges[edge];

    auto mu = matched_edge[u];
    auto mv = matched_edge[v];

    if (mu != NetworKit::none && edge_in_matching[mu]) remove_edge_from_matching(mu);
    if (mv != NetworKit::none && edge_in_matching[mv]) remove_edge_from_matching(mv);

    edge_in_matching[edge] = true;
    edges_in_matching ++;
    matched_vertex[u] = v;
    matched_vertex[v] = u;
    matched_edge[u] = matched_edge[v] = edge;
}

void GabowScalingMatching::remove_edge_from_matching(NetworKit::edgeid edge) {
    if (!edge_in_matching[edge]) {
        return;
    }

    auto [u, v] = graph_edges[edge];

    edge_in_matching[edge] = false;
    edges_in_matching --;
    if (matched_vertex[u] == v)
        matched_vertex[u] = matched_edge[u] = NetworKit::none;
    if (matched_vertex[v] == u)
        matched_vertex[v] = matched_edge[v] = NetworKit::none;
}   

void GabowScalingMatching::check_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v] = graph_edges[edge];
    
    if (edge_in_matching[edge]) {
        matched_vertex[u] = v;
        matched_vertex[v] = u;
    }
}

void GabowScalingMatching::create_trivial_blossoms(OldBlossom* T) {
    for (auto node : T->nodes) {
        Blossom* node_blossom = new Blossom {
            node, // base
            0, 0, 0, 0, 0, // z0, t_root, t_odd, t_even, Delta
            nullptr, {}, // parent, subblossoms
            free, no_edge, false, // label, backtrack_edge, visited
            {}, nullptr // list iterators, split list
        };
        trivial_blossom[node] = node_blossom;
        add_blossom(node_blossom, T);
    }

    for (auto sT : T->children) 
        create_trivial_blossoms(sT);
}

GabowScalingMatching::OldBlossom* 
GabowScalingMatching::turn_current_blossoms_into_old(const std::list<Blossom*>& subblossoms) {
    OldBlossom* T = new OldBlossom;
    T->size = 0;
    T->z = 0;
    T->heavy_path_parent = nullptr;
    T->heavy_child = nullptr;
    T->dissolved = false;

    for (auto b : subblossoms) {
        if (b->is_trivial()) {
            T->nodes.push_back(b->base);
            T->size ++;
        } else {
            std::list<Blossom*> subblossoms;
            for (auto [sb, e] : b->subblossoms) subblossoms.push_back(sb);
            auto sT = turn_current_blossoms_into_old(subblossoms); 
            sT->z = b->z0;
            T->size += sT->size;
            T->children.push_back(sT);
        }

        delete b;
    }

    return T;
}

void GabowScalingMatching::heavy_path_decomposition(OldBlossom* T, MaximumWeightMatching::intweight outer_dual) {
    // Divide the old blossom tree into heavy paths
    // Call path on the roots of these paths in postorder order

    for (auto child : T->children) {
        if (2 * child->size > T->size) {
            T->heavy_child = child;
            child->heavy_path_parent = T;
        }

        heavy_path_decomposition(child, outer_dual + T->z);

        if (child != T->heavy_child) {
            // The child is not on the path so it has been dissolved
            // Add the nodes and blossoms it contains
            T->nodes.splice(T->nodes.end(), std::move(child->nodes));
            T->shell_blossoms.splice(T->shell_blossoms.end(), std::move(child->shell_blossoms));

            delete child;
        }
    }

    T->children.clear();

    // A root of a heavy path has been found, dissolve it
    if (T->is_heavy_path_root()) {
        path(T, outer_dual);
    }
}

void GabowScalingMatching::path(OldBlossom* B, MaximumWeightMatching::intweight outer_dual) {
    outer_shells_dual = outer_dual;
    path_root = highest_undissolved = B;
    path_root->for_nodes([this] (NetworKit::node v) {
        vertex_path[v] = path_root;
    });
    
    while (true) {
        enumerate_shells();

        // Maximize the number of matched edges using tight edges
        augmentPaths();

        // Finish dissolving the path when all shells are dissolved or the only remaining shell
        // represents the entire graph and the matching is perfect, which ends the whole algorithm
        if (shells.size() == 0 || 
            (shells.size() == 1 && shells.front().second == old_root && matching_is_perfect())) {
            return;
        }

        // Count free nodes in each shell
        for (int i = 0; i < shells.size(); ++ i) {
            shells[i].first = free_nodes_in_shell(shells[i].second);
        }
        
        // Search shells in the order of decreasing number of free nodes
        std::sort(shells.begin(), shells.end(), std::greater<>());
        for (auto [free_nodes, shell] : shells) {
            if (free_nodes > 0 && !shell->dissolved && !shell->searched) {
                shell_search(shell);
            }
        }

        // Update the dual weights with distributions done after it was searched
        if (highest_undissolved != nullptr) {
            highest_undissolved->for_nodes([this] (NetworKit::node v) {                
                current_y[v] += distribution_so_far(current_shell[v]->shell_index + 1) - Delta[v];
            });
        }

        // Delete the dissolved shells except the path root
        for (auto [_, s] : shells) {
            if (s->dissolved && s != path_root) {
                delete s;
            }
        }
    }
}

int GabowScalingMatching::free_nodes_in_shell(OldBlossom* B) {
    int free_nodes = 0;
    for (auto v : B->nodes) {
        if (matched_vertex[v] == NetworKit::none)
            free_nodes ++;
    }
    return free_nodes;
}

void GabowScalingMatching::enumerate_shells() {
    shells.clear();
    current_shell_duals.clear();

    // Find all the undissolved shells on the current path
    for (auto shell = highest_undissolved; shell != nullptr; shell = shell->heavy_child) {
        shell->searched = false;
        for (auto v : shell->nodes) {
            current_shell[v] = shell;
            search_shell[v] = nullptr;
            Delta[v] = 0;
            y0[v] = current_y[v];
            union_find.reset(v, trivial_blossom[v]);
        };

        // Remember the blossom of each vertex before path augmentation
        for (auto B : shell->shell_blossoms) {
            B->for_nodes([this, B] (NetworKit::node v) {
                current_blossom[v] = B;
            });
            B->label = free;
        }

        shells.push_back({0, shell});
        current_shell_duals.push_back(shell->z);
    }

    // Clear the distributions
    shell_distribution.reset(shells.size());

    // Record the index of each shell, sum up duals of above shells
    for (int i = 0; i < shells.size(); ++ i)  {
        shells[i].second->shell_index = i;
        if (i > 0) current_shell_duals[i] += current_shell_duals[i-1];
    }
}

void GabowScalingMatching::change_blossom_base(Blossom* B, NetworKit::node new_base, Edge edge) {
    if (matched_edge[B->base] != NetworKit::none) {
        remove_edge_from_matching(matched_edge[B->base]);
    }

    swap_edges_on_even_path(B, B->base, new_base);
}

void GabowScalingMatching::swap_edges_on_even_path(Blossom* B, NetworKit::node u, NetworKit::node v) {
    if (B->is_trivial() || u == v) return;

    auto out_vertex = u == B->base ? v : u;
    B->base = out_vertex;
    auto out_blossoms = blossoms_containing(out_vertex, B);

    swap_edges_on_even_path(B, out_vertex, std::move(out_blossoms));
}

void GabowScalingMatching::swap_edges_on_even_path(
    Blossom* B, NetworKit::node out_vertex, std::list<Blossom*>&& out_blossoms) {

    if (B->is_trivial()) return;
    
    B->base = out_vertex;
    Blossom* out_blossom = out_blossoms.front();
    out_blossoms.pop_front();
    auto [pathA, pathB] = split_subblossoms(B->subblossoms, out_blossom);
    Blossom* base_blossom = std::get<Blossom*>(B->subblossoms.back());

    if (pathA.size() % 2 == 0) {
        for (auto [b, e] : pathA) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathB) check_edge_in_matching(e.id);

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, base_blossom->base, (pathA.begin()->second).u);

        for (auto iter = pathA.begin(); std::next(iter) != pathA.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(out_blossoms));
    } else {
        for (auto [b, e] : pathB) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathA) check_edge_in_matching(e.id);

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(out_blossoms));

        for (auto iter = pathB.begin(); std::next(iter) != pathB.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, base_blossom->base, pathB.back().second.v);
    }

    pathB.splice(pathB.end(), pathA);
    B->subblossoms = pathB; // pathB + pathA   
}

void GabowScalingMatching::augmentPaths() {
    std::vector<std::tuple<NetworKit::node, NetworKit::node, Edge>> edges;
    NetworKit::index counter = 0;

    // Contract the graph into blossoms
    path_root->for_nodes([this, &counter] (NetworKit::node v) {
        if (v == current_blossom[v]->base) {
            actual_to_contracted[v] = counter;
            counter ++;
        }
    });

    NetworKit::Graph T(counter, false, false, true);
    std::vector<NetworKit::node> initial_matching(counter, NetworKit::none);

    if (highest_undissolved == nullptr) return;

    highest_undissolved->for_nodes([this, &edges, &T, &initial_matching] (NetworKit::node v) {
        reducedGraph.forEdgesOf(v, 
            [this, &edges, &T, &initial_matching] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid e) {

                // Contract blossoms, only consider tight edges within the same shell
                if (vertex_path[u] != path_root || u < v || current_shell[v] != current_shell[u] ||
                    current_blossom[v] == current_blossom[u] || current_slack({u, v, e}) > 0) 
                    return;

                auto cv = actual_to_contracted[current_blossom[v]->base];
                auto cu = actual_to_contracted[current_blossom[u]->base];

                if (cu > cv) {
                    std::swap(cu, cv);
                    std::swap(u, v);
                }

                if (edge_in_matching[e]) {
                    initial_matching[cu] = cv;
                    initial_matching[cv] = cu;
                }

                T.addEdge(cu, cv);
                edges.push_back(std::make_tuple(cu, cv, Edge{u, v, e}));
        });
    });

    // Sort the edges to be able to find dthem
    std::sort(edges.begin(), edges.end());
    // Maximize the number of edges in the matching using a maximum cardinality matching algorithm
    MicaliVaziraniMatching cardinality_matching(T, initial_matching);
    cardinality_matching.run();
    auto matching = cardinality_matching.getMatching();

    // Look through all edges in a matching
    for (NetworKit::index cv = 0; cv < counter; ++ cv) {
        if (matching[cv] != NetworKit::none && cv < matching[cv]) {
            auto cu = matching[cv];

            // Find the shell in a list to identify it
            auto eit = std::prev(std::lower_bound(
                edges.begin(), edges.end(), 
                std::make_tuple(cv, cu, no_edge)
            ));
            auto v = std::get<Edge>(*eit).u;
            auto u = std::get<Edge>(*eit).v;
            auto e = std::get<Edge>(*eit).id;

            auto v_blossom = current_blossom[v];
            auto u_blossom = current_blossom[u];

            // Set the edge in matching 
            change_blossom_base(v_blossom, v, {v, u, e});
            change_blossom_base(u_blossom, u, {u, v, e});
            set_edge_in_matching(e);
        }
    }
}

MaximumWeightMatching::intweight GabowScalingMatching::current_slack(Edge edge) {
    auto [u, v, id] = edge;

    // Calculate the slack of the edge 
    // Include duals of old blossoms above the current path,
    // the blossom that contains both vertices and old blossoms that 
    return current_y[u] + current_y[v] - current_w[id] + outer_shells_dual +
        (current_blossom[u] == current_blossom[v] ? current_blossom[v]->z0 : 0) +
        current_shell_duals[std::min(current_shell[u]->shell_index, current_shell[v]->shell_index)];
}

void GabowScalingMatching::shell_search(OldBlossom* B) {
    active_shell = starting_shell = B;
    init_shell_search();

    // Process events until the search is done
    while (!search_done) {
        auto event = event_queue.getEvent();

        switch (event.type) {
            case Event::Type::grow:
                grow(event.args.grow.v, event.args.grow.e);
                break;
            case Event::Type::blossom:
                blossom(event.args.uv);
                break;
            case Event::Type::dissolve:
                dissolve(event.args.B);
                break;
            case Event::Type::dissolveShell:
                dissolveShell(event.args.S);
                break;
        }
    }
}

void GabowScalingMatching::init_shell_search() {
    search_done = false;
    event_queue.reset();
    t_active = 0;
    active_shell->searched = true;
    active_shell_initial_dual = current_shell_duals[active_shell->shell_index] - 
                                2 * distribution_so_far(active_shell->shell_index + 1);

    for (auto b : active_shell->shell_blossoms) {
        // Create a list for each blossom
        b->list = split_find_min.init(b->node_list(), b);
        b->Delta = 0;
        b->label = free;

        // Schedule a grow event for an exposed blossom
        if (is_exposed(b)) {
            event_queue.scheduleEvent(0, Event::make_grow(b->base, no_edge));
        }
    }

    // Set counters for starting nodes
    for (auto v : active_shell->nodes) {
        t_shell[v] = 0;
        search_shell[v] = starting_shell;
        Delta[v] = distribution_so_far(active_shell->shell_index + 1);
    };

    // If the shell is not the outermost one, schedule it's dissolution
    if (active_shell != old_root) {
        event_queue.scheduleEvent(active_shell->z / 2, Event::make_dissolveShell(active_shell));
        t_undissolvable = 1e9;
    } else {
        t_undissolvable = 0;
    }

    // Schedule a dissolution of the inner shell if it exists
    auto inner_shell = active_shell->heavy_child;
    if (inner_shell != nullptr) {
        t_inner = 0;
        event_queue.scheduleEvent(inner_shell->z / 2, Event::make_dissolveShell(inner_shell));
    }
}

void GabowScalingMatching::finish_shell_search() {
    search_done = true;

    for (auto v : active_shell->nodes) {
        // Save the current value of the dual variable
        current_y[v] = y(v);

        // Record the ditributions done to v when it was in the search 
        Delta[v] += (std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    }

    for (auto B : active_shell->shell_blossoms) {
        // Cleanup the lists
        delete_lists(B);

        // Save the current value of the dual variable
        B->z0 = z(B);
    }

    // Expand blossoms with 0 dual weight
    std::vector<Blossom*> to_expand;
    for (auto B : active_shell->shell_blossoms)
        if (B->z0 == 0)
            to_expand.push_back(B);
    for (auto B : to_expand) 
        expand_blossom(B, active_shell);

    // If the searched shell didn't dissolve and it's not the outermost shell,
    // updates the value of it's dual weight and register the distributed weights
    if (!active_shell->dissolved && active_shell != old_root) {
        // Calculate time the shell has been active
        auto distributed = event_queue.timeNow() - t_active;

        active_shell->z -= 2 * distributed;
        add_distribution(active_shell, distributed);
    }

    // Similarly update weights distributed by the innner shell
    auto inner_shell = active_shell->heavy_child;
    if (inner_shell != nullptr) {
        // Calculate the time the shell has been inner
        auto distributed = event_queue.timeNow() - t_inner;
        inner_shell->z -= 2 * distributed;
        add_distribution(inner_shell, distributed);
    }

    // Finally, mark the shell as searched
    active_shell->searched = true;
}

void GabowScalingMatching::delete_lists(Blossom* B) {
    if (B->list != nullptr) {
        split_find_min.deleteList(B->list);
        B->list = nullptr;

        // If the blossom has a list, none of it's children do
        return;
    }

    // Some children of the blossom might have associated lists
    for (auto [b, e] : B->subblossoms) {
        delete_lists(b);
    }
}

void GabowScalingMatching::grow(NetworKit::node v, Edge e) {
    // if (e != no_edge) union_find.addegde(e.u, e.v)

    // Check if the vertex is free
    auto B = split_find_min.list(v);
    if (B->label != free) 
        return;

    B->backtrack_edge = e;

    if (e != no_edge && slack(e) != 0) return;

    if (edge_in_matching[e.id] || e == no_edge) {
        // Even blossom found, record the time
        B->label = even;
        B->t_root = B->t_even = event_queue.timeNow();
        
        if (B->is_trivial()) {
            // Set the blossom of v to it's corresponding trivial blossom
            union_find.reset(v, B);

            // Schedule events associated with the even vertex
            schedule(v);
        } else {
            // Link the even vertices into one set with the value 
            B->for_nodes([this, B, v] (NetworKit::node u) {
                if (u == v) return;
                // union_find.addedge(v, u);
                union_find.link(v, u, B);
            });

            // Schedule events associated with even vertices
            B->for_nodes([this] (NetworKit::node u) {
                schedule(u);
            });
        }
    } else {
        // Odd blossom found, record the time
        B->label = odd;
        B->t_root = B->t_odd = event_queue.timeNow();

        // TODO call addedge on even path from v to B->base
        // Schedule events associated with the odd blossom
        schedule(B->base);
    }
}

void GabowScalingMatching::schedule(NetworKit::node u) {
    auto u_blossom = get_blossom(u);

    if (u_blossom->label == odd) {
        if (!u_blossom->is_trivial()) {
            // Schedule the expansion of an odd blosssom
            event_queue.scheduleEvent(
                event_queue.timeNow() + z(u_blossom) / 2, 
                Event::make_dissolve(u_blossom)
            );
        }

        // Continue the search through the vertex matched to u
        auto v = matched_vertex[u];
        auto e = matched_edge[u];
        auto v_blossom = split_find_min.list(v);

        // If v is free, continue growing the search tree, if not a blossom has been found
        if (v_blossom->label == free)
            event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, {u, v, e}));
        else 
            event_queue.scheduleEvent(event_queue.timeNow(), Event::make_blossom({v, u, e}));

    } else {
        // The vertex u is newly even
        // Scan the edges adjacent to it
        reducedGraph.forEdgesOf(u, [this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            // Ignore neighbours outside the shell, inside the same blossom or matched to u
            if (edge_in_matching[id] || vertex_path[v] != path_root || 
                search_shell[v] != starting_shell ||
                union_find.find(u) == union_find.find(v)) return;

            auto edge_slack = slack({u, v, id});
            auto v_blossom = split_find_min.list(v);

            if (v_blossom->label == even) {
                // Found an edge connecting two even blossoms
                // Once it becomes tight, a blossom will be formed
                // Schedule the blossom creation event based on the remaining slack
                event_queue.scheduleEvent(
                    event_queue.timeNow() + edge_slack / 2, 
                    Event::make_blossom({u, v, id})
                );
            } else if (v_blossom->label == odd) {
                // Check if the edge has a smaller slack than any other edge from an even vertex to v
                auto [key, min_edge] = split_find_min.currentKey(v);
                if (min_edge == no_edge || edge_slack < slack(min_edge)) {
                    // If v is in and odd blossom B, than 
                    // minimum slack(v, u) for even u is key(u) - (t_odd(B) - Delta(B))
                    split_find_min.decreaseKey(
                        v, edge_slack + (v_blossom->t_odd - v_blossom->Delta), {u, v, id});
                }
            } else if (v_blossom->label == free) {
                // If v is in a free blossom B, than 
                // minimum slack(v, u) for even u is key(u) - (t_now - Delta(B))
                auto min_slack = split_find_min.findMin(v_blossom->list).first - 
                                 (event_queue.timeNow() - v_blossom->Delta);
                
                split_find_min.decreaseKey(
                        v, edge_slack + (event_queue.timeNow() - v_blossom->Delta), {u, v, id});
                
                // Schedule a grow event if an edge with smaller slack has been found
                if (edge_slack < min_slack) {
                    event_queue.scheduleEvent(
                        event_queue.timeNow() + edge_slack, 
                        Event::make_grow(v, {u, v, id})
                    );
                }
            }
        });
    }
}

void GabowScalingMatching::dissolve(Blossom* B) {
    if (B->parent != nullptr || B->is_trivial()) return;

    // Dissolve the odd blossom 
    auto node = B->backtrack_edge.v;
    Blossom* inner_blossom = blossoms_containing(node, B).front();
    Blossom* base_blossom = B->subblossoms.back().first;
    auto [pathA, pathB] = split_subblossoms(B->subblossoms, inner_blossom);

    for (auto [b, e] : B->subblossoms) b->label = free;

    inner_blossom->label = odd;
    inner_blossom->backtrack_edge = B->backtrack_edge;

    if (pathB.size() % 2 == 0) {
        int parity = 0;
        for (auto iter = pathB.begin(); iter != pathB.end(); iter ++, parity ++) {
            iter->first->label = parity % 2 == 0 ? even : odd;
            iter->first->backtrack_edge = iter->second;
        }
    } else {
        int parity = 0;
        Blossom* prev = base_blossom;
        for (auto iter = pathA.begin(); iter != pathA.end(); iter ++, parity ++) {
            prev->label = parity % 2 == 0 ? odd : even;
            prev->backtrack_edge = reverse(iter->second);
            prev = iter->first;
        }
    }

    // Split the expanded blossom list of nodes to it's children
    for (auto it = B->subblossoms.begin(); std::next(it) != B->subblossoms.end(); it ++) {
        auto this_blossom = it->first;
        auto next_blossom = std::next(it)->first;
        auto [L1, L2] = split_find_min.split(this_blossom->base, this_blossom, next_blossom);
        this_blossom->list = L1;
        next_blossom->list = L2;
    }

    std::vector<NetworKit::node> to_schedule;
    for (auto [b, e] : B->subblossoms) {
        // Add b to the shell's list
        add_blossom(b, active_shell);
        
        // Update counters for b
        b->parent = nullptr;
        b->t_root = event_queue.timeNow();
        b->Delta = B->Delta + (event_queue.timeNow() - B->t_odd);

        if (b->label == odd) {
            b->t_odd  = event_queue.timeNow();

            // Schedule odd blossom expansion
            if (!b->is_trivial())
                event_queue.scheduleEvent(event_queue.timeNow() + b->z0 / 2, Event::make_dissolve(b));
        } else if (b->label == even) {
            b->t_even = event_queue.timeNow();

            // Link the even vertices
            b->for_nodes([this, b, &to_schedule] (NetworKit::node v) {
                if (v != b->base)
                    union_find.link(b->base, v, b);
                to_schedule.push_back(v);
            });
        } else {
            // The blossom is free. Schedule the grow event for the blossom by finding the edge 
            // with smallest slack that connects it to an even vertex.
            auto [key, edge] = split_find_min.findMin(b->list);
            if (edge != no_edge) {
                auto [u, v, id] = edge;
                auto slack = key - (event_queue.timeNow() - b->Delta);
                event_queue.scheduleEvent(event_queue.timeNow() + slack, Event::make_grow(v, edge));
            }
        }
    }

    // Schedule events for newly even vertices
    for (auto v : to_schedule) {
        schedule(v);
    }

    // Delete the blossom from the active shell's list
    delete_blossom(B, active_shell);
    delete B;
}

std::list<GabowScalingMatching::Blossom*> 
GabowScalingMatching::blossoms_containing(NetworKit::node u, Blossom* until) {
    std::list<Blossom*> res;
    Blossom* iter = trivial_blossom[u];
    while (iter != until) { 
        res.push_front(iter);
        iter = iter->parent;
    }
    return res;
}

std::pair<std::list<std::pair<GabowScalingMatching::Blossom*, GabowScalingMatching::Edge>>,
    std::list<std::pair<GabowScalingMatching::Blossom*, GabowScalingMatching::Edge>>>
GabowScalingMatching::split_subblossoms(
    std::list<std::pair<Blossom*, Edge>> sub_blossoms, Blossom* blossom) {
            
    auto iter = sub_blossoms.begin();
    while (iter->first != blossom) iter ++;
    return {{sub_blossoms.begin(), std::next(iter)}, {std::next(iter), sub_blossoms.end()}};
}

void GabowScalingMatching::blossom(Edge edge) {
    auto [u, v, id] = edge;

    // The edge is contained in a single even blossom
    if (union_find.find(u) == union_find.find(v)) return;

    // Backtrack from the blossoms
    backtrack(get_blossom(u), get_blossom(v), edge);
}

void GabowScalingMatching::dissolveShell(OldBlossom* S) {
    S->dissolved = true;
    S->z = 0;

    // Update dual weights distributed by the active shell
    auto distributed = event_queue.timeNow() - t_active;
    active_shell->z -= 2 * distributed;
    active_shell_initial_dual -= 2 * distributed;
    add_distribution(active_shell, distributed);
    t_active = event_queue.timeNow();
    
    if (S == active_shell) {
        // Find the shell that the active shell dissolves into
        auto outer_shell = S->heavy_path_parent;

        std::vector<Edge> edges_to_schedule;

        if (S == highest_undissolved) {
            // Outermost shell is dissolve, finish the search
            finish_shell_search();

            // Find the new highest undissolved shell
            while (highest_undissolved != nullptr && highest_undissolved->dissolved)
                highest_undissolved = highest_undissolved->heavy_child;
        } else {
            if (outer_shell->searched) {
                // The shell is dissolved into an already searched shell
                finish_shell_search();
            } else {
                // Update counters for new nodes
                for (auto v : outer_shell->nodes) {
                    t_shell[v] = event_queue.timeNow();
                    Delta[v] = distribution_so_far(S->shell_index);
                    search_shell[v] = starting_shell;
                }

                // Update counters and labels, create lists for new blossoms
                for (auto b : outer_shell->shell_blossoms) {
                    b->list = split_find_min.init(b->node_list(), b);
                    b->Delta = 0;
                    b->label = free;
                }

                for (auto v : outer_shell->nodes) {
                    // Schedule a grow event for a free vertex
                    if (matched_vertex[v] == NetworKit::none)
                        event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, no_edge));
                    
                    // Check for edges connecting to even vertices to schedule grow events
                    reducedGraph.forEdgesOf(v, 
                        [this, &edges_to_schedule] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
                            // Check if the neighbours is in the current shell
                            if (vertex_path[u] != path_root || search_shell[u] != starting_shell)
                                return;
                            if (is_even(u)) {
                                edges_to_schedule.push_back({u, v, id});
                            }
                        }
                    );
                }
            }
        }

        if (outer_shell == nullptr) return;

        // Add nodes and blossoms to the outer shell
        outer_shell->nodes.splice(outer_shell->nodes.end(), S->nodes);
        outer_shell->shell_blossoms.splice(outer_shell->shell_blossoms.end(), std::move(S->shell_blossoms));  

        // Remove the active shell from the path
        outer_shell->heavy_child = active_shell->heavy_child;
        if (outer_shell->heavy_child != nullptr)
            outer_shell->heavy_child->heavy_path_parent = outer_shell;

        // Change the active shell, update the counters and dual variables
        active_shell = outer_shell;
        t_active = event_queue.timeNow();
        active_shell_initial_dual = current_shell_duals[active_shell->shell_index] - 
                                    2 * distribution_so_far(active_shell->shell_index + 1);
        active_shell->searched = true;

        if (search_done) return;

        if (active_shell != old_root) {
            // If the new active shell is not the outermost shell, schedule it's dissolution
            event_queue.scheduleEvent(
                event_queue.timeNow() + active_shell->z / 2,
                Event::make_dissolveShell(active_shell)
            );
        } else {
            // The new active shell is the root, so it becomes undissolvable, record the time
            t_undissolvable = event_queue.timeNow();
        }

        // Schedule grow events for relevant edges
        for (auto e : edges_to_schedule)
            event_queue.scheduleEvent(
                event_queue.timeNow() + slack(e), 
                Event::make_grow(e.v, e)
            );
    } else {
        // The inner shell is getting dissolved
        
        // Record the distributions from the innner shell before this event
        add_distribution(S, event_queue.timeNow() - t_inner);
        t_inner = event_queue.timeNow();

        std::vector<Edge> edges_to_schedule;

        if (S->searched) {
            // The inner shell has already been searched, finish the search
            finish_shell_search();
        } else {
            // Update values for newly added nodes
            for (auto v : S->nodes) {
                t_shell[v] = event_queue.timeNow();
                Delta[v] = distribution_so_far(S->shell_index + 1);
                search_shell[v] = starting_shell;
            };

            // Update values for newly added blossoms
            for (auto b : S->shell_blossoms) {
                b->list = split_find_min.init(b->node_list(), b);
                b->Delta = 0;
                b->label = free;
            }

            // Schedule events for newly added no
            for (auto v : S->nodes) {
                // Add exposed vertices to the search tree
                if (matched_vertex[v] == NetworKit::none)
                    event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, no_edge));
                
                // Check for edges connecting to even vertices to schedule grow events
                reducedGraph.forEdgesOf(v, 
                    [this, &edges_to_schedule] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
                        // Check if the neighbours is in the current shell
                        if (vertex_path[u] != path_root || search_shell[u] != starting_shell)
                            return;
                        if (is_even(u))
                            edges_to_schedule.push_back({u, v, id});
                    }
                );
            }
        }
        
        // Add the nodes and blossom from the inner shell to the active one
        active_shell->nodes.splice(active_shell->nodes.end(), S->nodes);
        active_shell->shell_blossoms.splice(active_shell->shell_blossoms.end(), std::move(S->shell_blossoms));

        // Remove the inner shell from the path
        active_shell->heavy_child = S->heavy_child;

        if (active_shell->heavy_child != nullptr) {
            // There is a new inner shell, update counters
            active_shell->heavy_child->heavy_path_parent = active_shell;

            if (search_done) return;

            // Schedule the dissolution of the new inner shell
            event_queue.scheduleEvent(
                event_queue.timeNow() + active_shell->heavy_child->z / 2,
                Event::make_dissolveShell(active_shell->heavy_child)
            );
        }

        if (search_done) return;

        // Schedule grow events for relevant edges
        for (auto e : edges_to_schedule)
            event_queue.scheduleEvent(
                event_queue.timeNow() + slack(e), 
                Event::make_grow(e.v, e)
            );
    }
}

bool GabowScalingMatching::backtrack(Blossom* u, Blossom* v, Edge edge) {
    std::vector<BacktrackInfo> u_path;
    std::vector<BacktrackInfo> v_path;

    Blossom* u_iter = u;
    Blossom* v_iter = v;

    u->visited = true;
    v->visited = true;

    // Backtrack both paths simultaneously while marking visited blossoms
    // to avoid O(V) backtrack steps when new blossom is detected

    while (u_iter != nullptr || v_iter != nullptr) {
        if (backtrack_step(u_iter, u_path) || backtrack_step(v_iter, v_path)) {
            // The two path intersect - a new blossom is found

            // Reset visisted flags for all blossom on the two paths
            for (auto [blossom, e] : u_path) blossom->visited = false; 
            for (auto [blossom, e] : v_path) blossom->visited = false; 
            u->visited = false; 
            v->visited = false;

            create_new_blossom(u, v, edge, u_path, v_path);

            return false;
        }
    }

    // Reset visisted flags for all blossom on the two paths
    for (auto [blossom, e] : u_path) blossom->visited = false;
    for (auto [blossom, e] : v_path) blossom->visited = false;
    u->visited = false; 
    v->visited = false;

    // The two paths don't intersect - an augmenting path was found ending the stage
    finish_shell_search();

    return true;
}

bool GabowScalingMatching::backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path) {
    if (iter == nullptr) return false;

    auto [a, b, id] = iter->backtrack_edge;

    if (a == NetworKit::none) {
        iter = nullptr;
        return false;
    }

    iter = get_blossom(a);
    path.push_back({iter, {a, b, id}});

    if (iter->visited) return true;

    iter->visited = true;
    return false;
}

void GabowScalingMatching::create_new_blossom(
        Blossom* u, Blossom* v, Edge edge, 
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    
    // Exclude blossoms after the one where paths meet
    cut_path_at(u_path, v_path.size() > 0 ? v_path.back().blossom : v, v);
    cut_path_at(v_path, u_path.size() > 0 ? u_path.back().blossom : u, u);

    // Remove repeated appearances of blossom u and v
    if (v_path.size() > 0 && v_path.back().blossom == u) u_path.clear();
    if (u_path.size() > 0 && u_path.back().blossom == v) v_path.clear();

    // The base is at the end of the path or it's v if the path is empty
    Blossom* base = v_path.size() > 0 ? v_path.back().blossom : v;
    std::list<std::pair<Blossom*, Edge>> subblossoms;
    
    // Construct the subblossom list
    for (int i = v_path.size() - 1; i > 0; -- i)
        subblossoms.emplace_back(v_path[i-1].blossom, v_path[i].edge);
    if (v_path.size() > 0) subblossoms.emplace_back(v, v_path[0].edge);
    subblossoms.emplace_back(u, reverse(edge));
    for (int i = 0; i < u_path.size(); ++ i)
        subblossoms.emplace_back(u_path[i].blossom, reverse(u_path[i].edge));

    Blossom* new_blossom = new Blossom {
        base->base, // base
        0, event_queue.timeNow(), 0, event_queue.timeNow(), 0, // z0, t_root, t_odd, t_even, Delta
        nullptr, subblossoms, // parent, subblossoms
        even, base->backtrack_edge, false, // label, backtrack_edge, visited
        {}, nullptr // list iterators, split list
    };

    std::vector<NetworKit::node> to_schedule;
    for (auto [b, edge] : subblossoms) {
        b->z0 = z(b);
        b->parent = new_blossom;
        if (b->label == odd) {
            b->label = even;

            // Update the counters for b
            b->t_even = event_queue.timeNow();
            b->Delta = b->Delta + (event_queue.timeNow() - b->t_odd);

            b->for_nodes([this, b, new_blossom, &to_schedule] (NetworKit::node v) {
                // Schedule events for newly even vertices
                to_schedule.push_back(v);

                // Link all nodes to the base
                if (v != b->base)
                    union_find.link(b->base, v, new_blossom);
            });
        }

        // Link the base to the new base so that all vertices of the new blossom are in one set
        union_find.link(new_blossom->base, b->base, new_blossom);

        // Remove b from the active shell's list
        delete_blossom(b, active_shell);
    }

    // Add the new blossom to the active shells list
    add_blossom(new_blossom, active_shell);

    // Schedule events for all the vertices that have become even for the first time
    for (auto v : to_schedule) 
        schedule(v);
}

void GabowScalingMatching::cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut2) {
    int index = -1;
    for (unsigned i = 0; i < path.size(); ++ i) {
        if (path[i].blossom == cut || path[i].blossom == cut2) {
            index = i;
            break;
        }
    }
    if (index != -1) path.resize(index + 1);
}

MaximumWeightMatching::intweight GabowScalingMatching::y(NetworKit::node v) {
    // To calculate current dual variable start with the value before search
    // Add the distributions from old blossoms before and after it became active
    MaximumWeightMatching::intweight basic = y0[v] + Delta[v] + 
            std::max(0, std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    auto B = split_find_min.list(v);

    // Add the dual adjustments based on it's label
    if (B->label == odd) {
        return basic + B->Delta + event_queue.timeNow() - B->t_odd;
    } else if (B->label == even) {
        return basic + B->Delta - (event_queue.timeNow() - B->t_even);
    }
    return basic + B->Delta;
}

MaximumWeightMatching::intweight GabowScalingMatching::z(Blossom* B) {
    // Calculate the value of the blossom's dual variable
    if (B->label == odd) {
        return B->z0 - 2 * (event_queue.timeNow() - B->t_odd);
    } else if (B->label == even) {
        return B->z0 + 2 * (event_queue.timeNow() - B->t_even);
    } else {
        return B->z0;
    }
}

MaximumWeightMatching::intweight GabowScalingMatching::shell_z() {
    // Calculate the current weight of the active shell
    return active_shell == old_root ? 0 :
        (active_shell_initial_dual - 2 * (event_queue.timeNow() - t_active));
}

MaximumWeightMatching::intweight GabowScalingMatching::slack(Edge e) {
    auto [u, v, id] = e;
    auto u_B = get_blossom(u);
    auto v_B = get_blossom(v); 

    // Calculate the slack of the edge
    // Include the dual weights of all the shells that contain the edge
    return y(u) + y(v) - current_w[id] + (u_B == v_B ? z(u_B) : 0) + outer_shells_dual + shell_z();
}

bool GabowScalingMatching::is_exposed(Blossom *b) {
    return matched_vertex[b->base] == NetworKit::none;
}

bool GabowScalingMatching::is_even(NetworKit::node v) {
    return split_find_min.list(v)->label == even;
}

bool GabowScalingMatching::is_odd(NetworKit::node v) {
    return split_find_min.list(v)->label == odd;
}

void GabowScalingMatching::add_distribution(OldBlossom* S, MaximumWeightMatching::intweight distribution) {
    shell_distribution.add(S->shell_index, distribution);
}

MaximumWeightMatching::intweight GabowScalingMatching::distribution_so_far(int shell_index) {
    return shell_distribution.sum(shell_index - 1);
}

GabowScalingMatching::Blossom* GabowScalingMatching::get_blossom(NetworKit::node v) {
    auto B = split_find_min.list(v);

    if (B->label == even) {
        // If the vertex is even the blossom is maintained with union-find
        return union_find.find(B->base);
    } else {
        // For non-even vertices it's blossom corresponds to the list in split-findmin
        return B;
    }
}

bool GabowScalingMatching::OldBlossom::is_heavy_path_root() {
    return heavy_path_parent == nullptr;
}

void GabowScalingMatching::OldBlossom::for_nodes(const std::function<void(NetworKit::node)>& handle) {
    for (auto n : nodes)
        handle(n);
    for (auto b : children)
        b->for_nodes(handle);
    if (heavy_child != nullptr) 
        heavy_child->for_nodes(handle);
}

void GabowScalingMatching::OldBlossom::for_blossoms(const std::function<void(OldBlossom*)>& handle) {
    handle(this);
    for (auto b : children)
        b->for_blossoms(handle);
    if (heavy_child != nullptr) 
        heavy_child->for_blossoms(handle);
}

bool GabowScalingMatching::Blossom::is_trivial() {
    return subblossoms.size() == 0;
}

void GabowScalingMatching::Blossom::for_nodes(const std::function<void(NetworKit::node)>& handle) {
    if (is_trivial()) {
        handle(base);
    } else {
        for (auto [b, e] : subblossoms)
            b->for_nodes(handle);
    }
}

std::list<NetworKit::node> GabowScalingMatching::Blossom::node_list() {
    std::list<NetworKit::node> nodes;
    for_nodes([this, &nodes] (NetworKit::node v) {
        nodes.push_back(v);
    });
    return nodes;
}

bool GabowScalingMatching::matching_is_perfect() {
    return edges_in_matching == reducedGraph.numberOfNodes() / 2;
}

GabowScalingMatching::Edge GabowScalingMatching::reverse(const Edge& info) {
    auto [u, v, id] = info;
    return { v, u, id };
}

bool GabowScalingMatching::Edge::operator==(const GabowScalingMatching::Edge& other) const {
    return u == other.u && v == other.v && id == other.id;
}

bool GabowScalingMatching::Edge::operator!=(const GabowScalingMatching::Edge& other) const {
    return !(*this == other);
}

bool GabowScalingMatching::Edge::operator<(const GabowScalingMatching::Edge& other) const {
    return u  < other.u ? true  : 
           u  > other.u ? false :
           v  < other.v ? true  :
           v  > other.v ? false :
           id < other.id;
}

void GabowScalingMatching::Blossom::short_print() {
    if (is_trivial()) {
        std::cerr << base;
        return;
    }
    std::cerr << "{"; 
    for (auto [b, e] : subblossoms) {
        std::cerr << "(";
        b->short_print();
        std::cerr << ", (" << e.u << ", " << e.v << ")) ";
    }
    std::cerr << ":" << z0;
    if (parent == nullptr) {
        std::cerr << ":" << label_to_str(label);
    }
    std::cerr << "}";
}

void GabowScalingMatching::Blossom::nodes_print() {
    std::cerr << "("; 
    if (is_trivial()) {
        std::cerr << base;
    } else {
        for (auto sb : subblossoms) sb.first->nodes_print();
    }
    std::cerr << ")";
}

void GabowScalingMatching::OldBlossom::short_print() {
    std::cerr << "{ "; 
    for (auto v : nodes) {
        std::cerr << v << " ";
    }
    for (auto b : children) {
        b->short_print();
        std::cerr << " ";
    }
    if (heavy_child != nullptr) {
        std::cerr << " >";
        heavy_child->short_print();
        std::cerr << " ";
    }
    std::cerr << "} : " << z;
}

std::ostream& operator<<(std::ostream &out, const GabowScalingMatching::Edge& edge) {
    out << "(" << node_to_str(edge.u) << ", " << node_to_str(edge.v) << ")";
    return out;
}

std::ostream& operator<<(std::ostream &out, const GabowScalingMatching::Event& event) {
    switch (event.type) {
        case GabowScalingMatching::Event::Type::grow:
            out << "grow(" << event.args.grow.v << ", " << event.args.grow.e << ")";
            break;
        case GabowScalingMatching::Event::Type::blossom:
            out << "blossom(" << event.args.uv << ")";
            break;
        case GabowScalingMatching::Event::Type::dissolve:
            out << "dissolve(";
            event.args.B->nodes_print();
            out << ")";
            break;
        case GabowScalingMatching::Event::Type::dissolveShell:
            out << "dissolveShell(" << event.args.S << ")";
            break;
    }
    return out;
}

void GabowScalingMatching::print_backtrack(
        Blossom* u, Blossom* v, Edge edge, 
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    
    std::cerr << "u_path:\n";
    for (auto [b, e] : u_path) {
        std::cerr << "   >";
        b->short_print(); std::cerr << std::endl;
        std::cerr << "    (" << e.u << ", " << e.v << ")\n";
    }
    std::cerr << "u:  "; u->short_print(); std::cerr << std::endl;
    {
        auto [p, q, id] = edge;
        std::cerr << "    (" << p << ", " << q << ")\n";
    }
    std::cerr << "v:  "; v->short_print(); std::cerr << std::endl;
    std::cerr << "v_path:\n";
    for (auto [b, e] : v_path) {
        std::cerr << "   >";
        b->short_print(); std::cerr << std::endl;
        std::cerr << "    (" << e.u << ", " << e.v << ")\n";
    }
}

void GabowScalingMatching::print_search_state() {
    std::cerr << "shell dual: " << shell_z() << std::endl;
    std::cerr << "node: ";
    for (auto v : active_shell->nodes) std::cerr << " \t" << v;
    std::cerr << "\ny     ";
    for (auto v : active_shell->nodes) std::cerr << " \t" << y(v);
    std::cerr << "\nmatch ";
    for (auto v : active_shell->nodes) std::cerr << " \t" << node_to_str(matched_vertex[v]);
    std::cerr << "\nlabel ";
    for (auto v : active_shell->nodes) std::cerr << " \t" << label_to_str(get_blossom(v)->label);
    std::cerr << std::endl;
    std::cerr << "nontrivial blossoms:\n";
    for (auto b : active_shell->shell_blossoms) {
        if (b->is_trivial()) continue;
        b->short_print();
        std::cerr << " : " << z(b) << std::endl;
    }
    std::cerr << "edges:  ";
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        if (search_shell[u] == starting_shell && search_shell[v] == starting_shell) {
            std::cerr << "\t" << u << ", " << v << "";
        }
    });
    std::cerr << "\nslack:  ";
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        if (search_shell[u] == starting_shell && search_shell[v] == starting_shell) {
            auto s = slack({u, v, e});
            std::cerr << "\t" << s;
        }
    });
    std::cerr << std::endl;
    std::cerr << "shell_distributions: ";
    for (int i = 0; i < shells.size(); ++ i) 
        std::cerr << distribution_so_far(i + 1) << " ";
    std::cerr << std::endl;
    for (auto it = active_shell->shell_blossoms.begin(); it != active_shell->shell_blossoms.end(); it ++) {
        assert((*it)->shell_blossoms_it == it);
    }
}

void GabowScalingMatching::print_current_state() {
    std::cerr << "node: ";
    path_root->for_nodes([this] (NetworKit::node v) {std::cerr << "\t" << v; });
    std::cerr << "\ny:    ";
    path_root->for_nodes([this] (NetworKit::node v) { std::cerr << "\t" << current_y[v]; });
    std::cerr << "\nmatch ";
    path_root->for_nodes([this] (NetworKit::node v) { 
        std::cerr << "\t" << node_to_str(matched_vertex[v]); 
    });
    std::cerr << std::endl;
    std::cerr << "edges: ";
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        if (vertex_path[u] == path_root && vertex_path[v] == path_root) {
            std::cerr << "\t" << u << ", " << v;
        }
    });
    std::cerr << "\nslack: ";
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        if (vertex_path[u] == path_root && vertex_path[v] == path_root) {
            auto s = current_slack({u, v, e});
            std::cerr << "\t" << s;
        }
    });
    std::cerr << std::endl;
    for (auto [_, shell] : shells) {
        for (auto it = shell->shell_blossoms.begin(); it != shell->shell_blossoms.end(); it ++) {
            assert((*it)->shell_blossoms_it == it);
        }
    }
    int odd_free = 0;
    int even_free = 0;
    path_root->for_nodes([this, &odd_free, &even_free] (NetworKit::node v) { 
        if (matched_vertex[v] != NetworKit::none) return;
        if (current_y[v] % 2 == 0)
            ++ even_free;
        else 
            ++ odd_free;
    });
    std::cerr << odd_free << " " << even_free << std::endl;
    assert(odd_free == 0 || even_free == 0);
}

void GabowScalingMatching::print_vertex_state(NetworKit::node v) {
    if (search_shell[v] != starting_shell) return;

    MaximumWeightMatching::intweight basic = y0[v] + Delta[v] + 
            std::max(0, std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    auto B = split_find_min.list(v);
    
    MaximumWeightMatching::intweight y;
    
    if (B->label == odd) {
        y = basic + B->Delta + event_queue.timeNow() - B->t_odd;
    } else if (B->label == even) {
        y = basic + B->Delta - (event_queue.timeNow() - B->t_even);
    } else {
        y = basic + B->Delta;
    }

    std::cerr << "Y(" << v << ") = " << y << ":\n";
    std::cerr << "y0: " << y0[v] << std::endl;
    std::cerr << "Delta: " << Delta[v] << std::endl;
    std::cerr << "t_shell: " << t_shell[v] << std::endl;
    std::cerr << "basic: " << basic << std::endl;
    std::cerr << "t_undissolvable: " << t_undissolvable << std::endl;
    std::cerr << "label: " << label_to_str(B->label) << std::endl;
    std::cerr << "B->Delta: " << B->Delta << std::endl;
    std::cerr << "B->t_odd: " << B->t_odd << std::endl;
    std::cerr << "B->t_even: " << B->t_even << std::endl;
}

std::string GabowScalingMatching::label_to_str(Label label) {
    switch (label)
    {
    case free:
        return "-";
        break;
    case even:
        return "S";
        break;
    case odd:
        return "T";
        break;
    default:
        break;
    }
}

void GabowScalingMatching::check_consistency() {
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        if (edge_in_matching[id]) {
            if ((matched_vertex[u] != v) || (matched_vertex[v] != u) ||
                (matched_edge[u] != id) || (matched_edge[v] != id)) {
                std::cerr << "Wrong values for matched edge (" << u << ", " << v << ")\n"; 
                assert(false);
            }
        } else {
            if ((matched_vertex[u] == v) || (matched_vertex[v] == u) ||
                (matched_edge[u] == id) || (matched_edge[v] == id)) {
                std::cerr << "Wrong values for unmatched edge (" << u << ", " << v << ")\n"; 
                assert(false);
            }
        }
    });
    path_root->for_blossoms([this] (OldBlossom* S) {
        for (auto B : S->shell_blossoms) {
            check_consistency(B);
        }
    });
}

void GabowScalingMatching::check_consistency(Blossom* B) {
    bool in = false;
    for (auto [b, e] : B->subblossoms) {
        check_consistency(b);
        if (edge_in_matching[e.id] != in) {
            std::cerr << "Wrong edge matching status " << e << " in blossom ";
            B->short_print();
            std::cerr << std::endl;
            assert(false);
        }
        in = !in;
    }
}

} /* namespace Koala */