#include <matching/MaximumMatching.hpp>

namespace Koala {

std::string node_to_str(NetworKit::node v) {
    return v == NetworKit::none ? "-" : std::to_string(v);
}


NetworKit::Graph reduce_to_MWPM(const NetworKit::Graph& graph) {
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
    MaximumMatching(graph), 
    reducedGraph(reduce_to_MWPM(graph)), 
    current_y(reducedGraph.upperNodeIdBound()),
    trivial_blossom(reducedGraph.upperNodeIdBound(), nullptr),
    union_find(reducedGraph.upperNodeIdBound()),
    matched_vertex(reducedGraph.upperNodeIdBound()),
    matched_edge(reducedGraph.upperNodeIdBound()),
    split_find_min(reducedGraph.upperNodeIdBound(), 1000000000, no_edge),
    y0(reducedGraph.upperNodeIdBound()),
    Delta(reducedGraph.upperNodeIdBound()),
    t_shell(reducedGraph.upperNodeIdBound()),
    current_blossom(reducedGraph.upperNodeIdBound(), nullptr),
    current_shell(reducedGraph.upperNodeIdBound()),
    search_shell(reducedGraph.upperNodeIdBound()),
    graph_edges(reducedGraph.upperEdgeIdBound()),
    actual_to_contracted(reducedGraph.upperEdgeIdBound()),
    edge_in_matching(reducedGraph.upperEdgeIdBound()),
    vertex_path(reducedGraph.upperNodeIdBound()) { 
        reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
            graph_edges[e] = {u, v};
        });
    }

void GabowScalingMatching::run() {
    std::vector<MaximumMatching::intedgeweight> w(reducedGraph.upperEdgeIdBound());
    
    reducedGraph.forEdges(
        [&w] (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight ew, NetworKit::edgeid e) {
            w[e] = 2 * static_cast<MaximumMatching::intedgeweight>(ew);
        }
    );

    #if DEBUG_LOGGING
    std::cerr << "ORIGINAL:\n";
    std::cerr << "NODES:\n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << " : ";
        graph.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u) {
            std::cerr << u << " ";
        });
        std::cerr << std::endl;
    });
    std::cerr << "EDGES:\n";
    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        std::cerr << u << " " << v << "\n";
    });

    std::cerr << "REDUCED:\n";
    std::cerr << "NODES:\n";
    reducedGraph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << " : ";
        reducedGraph.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u) {
            std::cerr << u << " ";
        });
        std::cerr << std::endl;
    });
    std::cerr << "EDGES:\n";
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        std::cerr << e << " " << "(" << u << ", " << v << ")\n";
    });
    #endif

    auto [M, _, __] = scale(w);

    graph.forNodes([this, &M] (NetworKit::node v) {
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

    std::vector<int> w1(w.size());
    for (int i = 0; i < w.size(); ++ i) w1[i] = 2 * (w[i] / 4);

    auto [_, y1, T] = scale(w1);

    current_w = w;
    T->for_blossoms([this] (OldBlossom* B) { B->z = 2 * B->z; });
    reducedGraph.forNodes([this, &y1] (NetworKit::node v) {
        current_y[v] = 2 * y1[v] + 1;
        matched_vertex[v] = NetworKit::none;
        matched_edge[v] = NetworKit::none;
    });
    reducedGraph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        edge_in_matching[e] = false;
    });
    edges_in_matching = 0;
    old_root = T;

    #if DEBUG_LOGGING
    std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cerr << "SCALE\n";
    std::cerr << "W:\n";
    // for (int i = 0; i < current_w.size(); ++ i)
    //     std::cerr << i << " : (" << graph_edges[i].first << ", " 
    //               << graph_edges[i].second << ") : " << current_w[i] << std::endl;
    for (int i = 0; i < current_w.size(); ++ i) {
        auto [u, v] = graph_edges[i];
        if (u < graph.numberOfNodes() && v < graph.numberOfNodes()) {
            std::cerr << u << " " << v << "\t" << current_w[i] << std::endl;
        }
    }
    std::cerr << "Y:\n";
    for (int i = 0; i < current_y.size(); ++ i)
        std::cerr << i << " : " << current_y[i] << std::endl;
    std::cerr << "OLD BLOSSOM:\n";
    T->short_print();
    std::cerr << std::endl;
    #endif
    
    match(T);

    #if DEBUG_LOGGING
    T->shell_blossoms.sort([] (Blossom* a, Blossom* b) {
        return a->base < b->base;
    });
    std::cerr << "FINAL BLOSSOMS:\n";
    for (auto B : T->shell_blossoms) {
        B->short_print();
        std::cerr << std::endl;
    }
    #endif

    auto new_T = turn_current_blossoms_into_old(T->shell_blossoms);

    return std::make_tuple(matched_vertex, current_y, new_T);
}

void GabowScalingMatching::match(OldBlossom* T) {
    create_trivial_blossoms(T);
    
    heavy_path_decomposition(T, 0);
}

void GabowScalingMatching::delete_blossom(Blossom* B, OldBlossom* S) {
    #if DEBUG_LOGGING
    // std::cerr << "DELETE BLOSSOM " << B << " FROM " << S << std::endl;
    // for (auto b : S->shell_blossoms) std::cerr << b << " ";
    // std::cerr << std::endl;
    #endif

    S->shell_blossoms.erase(B->shell_blossoms_it);
}

void GabowScalingMatching::add_blossom(Blossom* B, OldBlossom* S) {
    #if DEBUG_LOGGING
    // std::cerr << "ADD BLOSSOM " << B << " TO " << S << std::endl;
    // for (auto b : S->shell_blossoms) std::cerr << b << " ";
    // std::cerr << std::endl;
    #endif

    S->shell_blossoms.push_back(B);
    B->shell_blossoms_it = std::prev(S->shell_blossoms.end());
}

void GabowScalingMatching::expand_blossom(Blossom* B, OldBlossom* S) {
    if (B->is_trivial()) return;
    
    for (auto [b, e] : B->subblossoms) {
        b->parent = nullptr;
        add_blossom(b, S);
    }

    delete_blossom(B, S);
    delete B;
}

void GabowScalingMatching::swap_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v] = graph_edges[edge];
    
    if (edge_in_matching[edge]) {
        #if DEBUG_LOGGING
        std::cerr << "Remove edge (" << u << ", " << v << ")" << std::endl;
        #endif

        edge_in_matching[edge] = false;
        edges_in_matching --;
        if (matched_vertex[u] == v)
            matched_vertex[u] = matched_edge[u] = NetworKit::none;
        if (matched_vertex[v] == u)
            matched_vertex[v] = matched_edge[v] = NetworKit::none;
    } else {
        #if DEBUG_LOGGING
        std::cerr << "Add edge (" << u << ", " << v << ")" << std::endl;
        #endif

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

    if (mu != NetworKit::none && edge_in_matching[mu]) swap_edge_in_matching(mu);
    if (mv != NetworKit::none && edge_in_matching[mv]) swap_edge_in_matching(mv);

    #if DEBUG_LOGGING
    std::cerr << "Set edge (" << u << ", " << v << ")\n";
    #endif

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

    #if DEBUG_LOGGING
    std::cerr << "Remove edge (" << u << ", " << v << ")\n";
    #endif

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
        #if DEBUG_LOGGING
        std::cerr << "Ensure edge (" << u << ", " << v << ")" << std::endl;
        #endif

        matched_vertex[u] = v;
        matched_vertex[v] = u;
    }
}

void GabowScalingMatching::create_trivial_blossoms(OldBlossom* T) {
    for (auto node : T->nodes) {
        Blossom* node_blossom = new Blossom {
            node, node, // base, initial_base
            0, 0, 0, 0, 0, // z0, t_root, t_odd, t_even, Delta
            nullptr, {}, // parent, subblossoms
            free, no_edge, false, // label, backtrack_edge, visited
            {}, nullptr // list iterators, split list
        };
        trivial_blossom[node] = node_blossom;
        add_blossom(node_blossom, T);
    }

    for (auto sT : T->subblossoms) 
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
            T->subblossoms.push_back(sT);
        }

        // delete b;
    }

    return T;
}

void GabowScalingMatching::heavy_path_decomposition(OldBlossom* T, MaximumMatching::intedgeweight outer_dual) {
    for (auto child : T->subblossoms) {
        if (2 * child->size > T->size) {
            T->heavy_child = child;
            child->heavy_path_parent = T;
        }

        heavy_path_decomposition(child, outer_dual + T->z);

        if (child != T->heavy_child) {
            T->nodes.splice(T->nodes.end(), std::move(child->nodes));
            T->shell_blossoms.splice(T->shell_blossoms.end(), std::move(child->shell_blossoms));

            delete child;
        }
    }

    T->subblossoms.clear();

    if (T->is_heavy_path_root()) {
        path(T, outer_dual);
    }
}

void GabowScalingMatching::path(OldBlossom* B, MaximumMatching::intedgeweight outer_dual) {
    #if DEBUG_LOGGING
    std::cerr << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n";
    std::cerr << "PATH ";
    B->short_print();
    std::cerr << " WITH " << outer_dual << std::endl;
    #endif

    outer_blossom_dual = outer_dual;
    path_root = highest_undissolved = B;
    path_root->for_nodes([this] (NetworKit::node v) {
        vertex_path[v] = path_root;
    });
    
    while (true) {
        #if DEBUG_LOGGING
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        #endif

        enumerate_shells();

        #if DEBUG_LOGGING
        print_current_state();
        #endif

        augmentPaths();

        for (int i = 0; i < shells.size(); ++ i) {
            shells[i].first = free_nodes_in_shell(shells[i].second);
        }
        std::sort(shells.begin(), shells.end(), std::greater<>());

        #if DEBUG_LOGGING
        for (auto [f, s] : shells) {
            std::cerr << s << " HAS " << f << " FREE NODES\n";
        }
        #endif

        if (shells.size() == 0 || (shells.size() == 1 && shells.front().second == old_root && matching_is_perfect())) {
            return;
        }
        
        for (auto [f, s] : shells) {
            if (!s->dissolved && f > 0 && !s->searched) {
                shell_search(s);
            }
        }

        #if DEBUG_LOGGING
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        std::cerr << "PATH ITERATION COMPLETE\n";
        #endif

        if (highest_undissolved != nullptr) {
            highest_undissolved->for_nodes([this] (NetworKit::node v) {
                #if DEBUG_LOGGING
                std::cerr << "y[" << v << "] += " << distribution_so_far(current_shell[v]->shell_index + 1) << " - " << Delta[v] << std::endl;
                #endif
                
                current_y[v] += distribution_so_far(current_shell[v]->shell_index + 1) - Delta[v];
            });
        }

        for (auto [_, s] : shells) {
            if (s->dissolved && s != path_root) {
                #if DEBUG_LOGGING
                std::cerr << "DELETING OLD BLOSSOM " << s << std::endl;
                #endif

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
    #if DEBUG_LOGGING
    std::cerr << "CURRENT SHELLS\n";
    if (path_root->dissolved) {
        assert(highest_undissolved == path_root->heavy_child);
    } else {
        assert(highest_undissolved == path_root);
    }
    #endif

    shells.clear();
    current_shell_duals.clear();
    // shell_distribution.clear();

    for (auto shell = highest_undissolved; shell != nullptr; shell = shell->heavy_child) {
        #if DEBUG_LOGGING
        std::cerr << "SHELL " << shell << " : " << shell->z << " : ";
        for (auto v : shell->nodes) std::cerr << v << " ";
        std::cerr << "\nBLOSSOMS:\n";
        for (auto B : shell->shell_blossoms) {
            B->short_print();
            if (B->is_trivial())
                std::cerr << std::endl;
            else std::cerr << " : " << B->z0 << std::endl;
        }
        shell->nodes.sort();
        #endif

        shell->searched = false;
        for (auto v : shell->nodes) {
            current_shell[v] = shell;
            search_shell[v] = nullptr;
            Delta[v] = 0;
            y0[v] = current_y[v];
            union_find.reset(v, trivial_blossom[v]);
        };
        for (auto B : shell->shell_blossoms) {
            B->for_nodes([this, B] (NetworKit::node v) {
                current_blossom[v] = B;
            });
            B->label = free;
        }

        if (shell->dissolved) continue;

        shells.push_back({0, shell});
        current_shell_duals.push_back(shell->z);
        // shell_distribution.push_back(0);
    }

    shell_distribution.reset(shells.size());

    for (int i = 0; i < shells.size(); ++ i)  {
        shells[i].second->shell_index = i;
        if (i > 0) current_shell_duals[i] += current_shell_duals[i-1];
    }

    // #if DEBUG_LOGGING
    // std::cerr << "CURRENT BLOSSOMS:\n";
    // path_root->for_nodes([this] (NetworKit::node v) {
    //     std::cerr << v << " : ";
    //     current_blossom[v]->short_print();
    //     std::cerr << std::endl;
    // });
    // #endif
}

void GabowScalingMatching::change_blossom_base(Blossom* B, NetworKit::node new_base, EdgeInfo edge) {
    // #if DEBUG_LOGGING
    // std::cerr << "CHANGE BASE TO " << new_base << " OF ";
    // B->short_print(); 
    // std::cerr << " FROM " << B->base << std::endl;
    // #endif

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

    // #if DEBUG_LOGGING
    // std::cerr << "Swap edges on even path " << B->base << " ~ " << out_vertex << " inside ";
    // B->nodes_print(); std::cerr << std::endl;
    // for (auto b : out_blossoms) {
    //     b->nodes_print();
    //     std::cerr << " ";
    // }
    // std::cerr << std::endl;
    // #endif

    if (B->is_trivial()) return;
    
    B->base = out_vertex;
    Blossom* out_blossom = out_blossoms.front();
    out_blossoms.pop_front();
    auto [pathA, pathB] = split_subblossoms(B->subblossoms, out_blossom);
    Blossom* base_blossom = std::get<Blossom*>(B->subblossoms.back());

    #if DEBUG_LOGGING
    std::cerr << "Augment to " << out_vertex << " inside "; out_blossom->nodes_print(); 
    std::cerr << " in blossom "; B->nodes_print(); std::cerr<<std::endl;
    #endif

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
    #if DEBUG_LOGGING
    std::cerr << "AUGMENT PATHS\n";
    #endif

    std::vector<std::tuple<NetworKit::node, NetworKit::node, EdgeInfo>> edges;
    NetworKit::index counter = 0;

    path_root->for_nodes([this, &counter] (NetworKit::node v) {
        if (v == current_blossom[v]->base) {
            // #if DEBUG_LOGGING
            // std::cerr << "CONTRACT ";
            // current_blossom[v]->short_print();
            // std::cerr << " INTO " << counter << std::endl;
            // #endif

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
                if (vertex_path[u] != path_root || u < v ||
                      current_shell[v] !=   current_shell[u] ||
                    current_blossom[v] == current_blossom[u] ||
                    current_slack({u, v, e}) > 0) return;

                auto cv = actual_to_contracted[current_blossom[v]->base];
                auto cu = actual_to_contracted[current_blossom[u]->base];

                if (edge_in_matching[e]) {
                    initial_matching[cu] = cv;
                    initial_matching[cv] = cu;
                }

                T.addEdge(cu, cv);
                if (cu < cv)
                    edges.push_back(std::make_tuple(cu, cv, EdgeInfo{u, v, e}));
                else
                    edges.push_back(std::make_tuple(cv, cu, EdgeInfo{v, u, e}));
        });
    });

    std::sort(edges.begin(), edges.end());
    MicaliVaziraniMatching cardinality_matching(T, initial_matching);
    cardinality_matching.run();
    auto matching = cardinality_matching.getMatching();

    #if DEBUG_LOGGING
    // std::cerr << "Contracted graph:\n";
    // T.forNodes([this, T] (NetworKit::node v) {
    //     std::cerr << v << " : ";
    //     T.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u) {
    //         std::cerr << u << " ";
    //     });
    //     std::cerr << std::endl;
    // });
    // std::cerr << "Initial matching:\n";
    // T.forNodes([this, &initial_matching] (NetworKit::node v) {
    //     std::cerr << v << " " << node_to_str(initial_matching[v]) << std::endl;
    // });
    std::cerr << "Matching:\n";
    T.forNodes([this, T, matching] (NetworKit::node v) {
        std::cerr << v << " " << node_to_str(matching.at(v)) << std::endl;
    });
    std::cerr << "Translated edges:\n";
    for (auto [cu, cv, e] : edges) {
        std::cerr << "(" << cu << ", " << cv << ") ~ " << e << "\n";
    }
    std::cerr << "Current blossoms:\n";
    for (NetworKit::index v = 0; v < reducedGraph.numberOfNodes(); v ++) {
        std::cerr << v << " : ";
        current_blossom[v]->short_print();
        std::cerr << std::endl;
    }
    #endif

    for (NetworKit::index cv = 0; cv < counter; ++ cv) {
        if (matching[cv] != NetworKit::none && cv < matching[cv]) {
            auto cu = matching[cv];

            auto eit = std::prev(std::lower_bound(
                edges.begin(), edges.end(), 
                std::make_tuple(cv, cu, no_edge)
            ));
            auto v = std::get<EdgeInfo>(*eit).u;
            auto u = std::get<EdgeInfo>(*eit).v;
            auto e = std::get<EdgeInfo>(*eit).id;

            auto v_blossom = current_blossom[v];
            auto u_blossom = current_blossom[u];

            #if DEBUG_LOGGING
            std::cerr << "SETTING EDGE (" << cv << ", " << cu
                      << ") ~ (" << v << ", " << u << ") BETWEEN BLOSSOMS ";
            v_blossom->short_print();
            std::cerr << " AND ";
            u_blossom->short_print();
            std::cerr << std::endl;
            #endif

            change_blossom_base(v_blossom, v, {v, u, e});
            change_blossom_base(u_blossom, u, {u, v, e});
            set_edge_in_matching(e);
        }
    }

    #if DEBUG_LOGGING
    std::cerr << "Matched edges " << edges_in_matching << std::endl;
    check_consistency();
    #endif
}

MaximumMatching::intedgeweight GabowScalingMatching::current_slack(EdgeInfo edge) {
    auto [u, v, id] = edge;

    // #if DEBUG_LOGGING
    //     std::cerr << "SLACK " << u << ", " << v << std::endl;
    //     std::cerr << current_y[u] << " " << current_y[v] << " " << current_w[id] << std::endl;
    //     std::cerr << current_blossom[u] << " " << current_blossom[v] << " " << current_blossom[v]->z0 << std::endl;
    //     std::cerr << current_shell_duals[std::min(current_shell[u]->shell_index, current_shell[v]->shell_index)] << std::endl;
    // #endif

    return current_y[u] + current_y[v] - current_w[id] + outer_blossom_dual +
        (current_blossom[u] == current_blossom[v] ? current_blossom[v]->z0 : 0) +
        current_shell_duals[std::min(current_shell[u]->shell_index, current_shell[v]->shell_index)];
}

void GabowScalingMatching::shell_search(OldBlossom* B) {
    #if DEBUG_LOGGING
    std::cerr << "==========================================================================================================\n";
    std::cerr << "SEARCH ";
    B->short_print();
    std::cerr << std::endl;
    #endif

    active_shell = starting_shell = B;
    init_shell_search();

    while (!search_done) {
        auto event = event_queue.getEvent();

        #if DEBUG_LOGGING
        std::cerr << "--------------------------------------------------------------------------\n";
        std::cerr << "EVENT " << event << " AT " << event_queue.timeNow() << std::endl;
        print_search_state();
        #endif

        switch (event.type) {
            case Event::Type::grow:
                grow(event.args.v, event.args.e);
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
    active_shell->searched = true;
    active_shell_initial_dual = current_shell_duals[active_shell->shell_index] - 
                                2 * distribution_so_far(active_shell->shell_index + 1);

    for (auto b : active_shell->shell_blossoms) {
        b->list = split_find_min.init(b->node_list(), b);
        if (is_exposed(b)) {
            event_queue.scheduleEvent(0, Event::make_grow(b->base, no_edge));
        }
        b->Delta = 0;
        b->label = free;

        // #if DEBUG_LOGGING
        // std::cerr << "Making list for blossom ";
        // b->short_print();
        // std::cerr << " -> " << b->list << std::endl;
        // #endif
    }

    for (auto v : active_shell->nodes) {
        t_shell[v] = 0;
        search_shell[v] = starting_shell;
        Delta[v] = distribution_so_far(active_shell->shell_index + 1);
    };

    active_shell->t_active = 0;
    if (active_shell != old_root) {
        event_queue.scheduleEvent(active_shell->z / 2, Event::make_dissolveShell(active_shell));
        t_undissolvable = 1e9;
    } else {
        t_undissolvable = 0;
    }

    auto inner_shell = active_shell->heavy_child;
    if (inner_shell != nullptr) {
        inner_shell->t_inner = 0;
        event_queue.scheduleEvent(inner_shell->z / 2, Event::make_dissolveShell(inner_shell));
    }
}

void GabowScalingMatching::finish_shell_search() {
    #if DEBUG_LOGGING
    std::cerr << "FINISH SHELL SEARCH\n";
    print_search_state();
    #endif

    search_done = true;

    for (auto B : active_shell->shell_blossoms) {
        delete_lists(B);
        B->z0 = z(B);
    }

    for (auto v : active_shell->nodes) {
        current_y[v] = y(v);
        Delta[v] += (std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    }

    std::vector<Blossom*> to_expand;
    for (auto B : active_shell->shell_blossoms)
        if (B->z0 == 0)
            to_expand.push_back(B);
    for (auto B : to_expand) expand_blossom(B, active_shell);

    if (!active_shell->dissolved && active_shell != old_root) {
        auto distributed = event_queue.timeNow() - active_shell->t_active;
        active_shell->z -= 2 * distributed;
        add_distribution(active_shell, distributed);

        #if DEBUG_LOGGING
        std::cerr << "z(" << active_shell << ") -= " << 2 * distributed << " = " << active_shell->z << std::endl;
        #endif
    }
    auto inner_blossom = active_shell->heavy_child;
    if (inner_blossom != nullptr) {
        auto distributed = event_queue.timeNow() - inner_blossom->t_inner;
        inner_blossom->z -= 2 * distributed;
        add_distribution(inner_blossom, distributed);

        #if DEBUG_LOGGING
        std::cerr << "z(" << inner_blossom << ") -= " << 2 * distributed << " = " << inner_blossom->z << std::endl;
        #endif
    }

    active_shell->searched = true;

    #if DEBUG_LOGGING
    std::cerr << "shell_distributions: ";
    for (int i = 0; i < shells.size(); ++ i) std::cerr << distribution_so_far(i + 1) << " ";
    std::cerr << std::endl;
    #endif
}

void GabowScalingMatching::delete_lists(Blossom* B) {
    // #if DEBUG_LOGGING
    // std::cerr << "DELETE LISTS ";
    // B->short_print();
    // std::cerr << std::endl;
    // #endif

    if (B->list != nullptr) {
        split_find_min.deleteList(B->list);
        B->list = nullptr;
        return;
    }

    for (auto [b, e] : B->subblossoms) {
        delete_lists(b);
    }
}

void GabowScalingMatching::grow(NetworKit::node v, EdgeInfo e) {
    // if (e != no_edge) union_find.addegde(e.u, e.v)

    auto B = split_find_min.list(v)->id;
    if (B->label != free) return;
    B->backtrack_edge = e;

    if (e != no_edge && slack(e) != 0) return; // TODO cancel before?

    #if DEBUG_LOGGING
    assert(e == no_edge || slack(e) == 0);
    #endif

    if (edge_in_matching[e.id] || e == no_edge) {
        B->label = even;
        B->t_root = B->t_even = event_queue.timeNow();
        
        if (B->is_trivial()) {
            union_find.reset(v, B);
            schedule(v);
        } else {
            B->for_nodes([this, B, v] (NetworKit::node u) {
                if (u == v) return;
                // union_find.addedge(v, u);
                union_find.link(v, u, B);
            });

            B->for_nodes([this] (NetworKit::node u) {
                schedule(u);
            });
        }
    } else {
        B->label = odd;
        B->t_root = B->t_odd = event_queue.timeNow();

        // TODO call addedge on even path from v to B->base

        schedule(B->base);
    }
}

void GabowScalingMatching::schedule(NetworKit::node u) {
    auto u_blossom = get_blossom(u);

    #if DEBUG_LOGGING
    std::cerr << "SCHEDULE FOR " << u << " INSIDE ";
    u_blossom->short_print();
    std::cerr << std::endl;
    #endif

    if (u_blossom->label == odd) {
        if (!u_blossom->is_trivial()) {
            assert(z(u_blossom) % 2 == 0);
            event_queue.scheduleEvent(
                event_queue.timeNow() + z(u_blossom) / 2, 
                Event::make_dissolve(u_blossom)
            );
        }

        auto v = matched_vertex[u];
        auto e = matched_edge[u];
        auto v_blossom = split_find_min.list(v)->id;

        assert(slack({u, v, e}) == 0);

        #if DEBUG_LOGGING
        std::cerr << u << " " << e << " ";
        v_blossom->short_print();
        std::cerr << std::endl;
        #endif

        if (v_blossom->label == free)
            event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, {u, v, e}));
        else 
            event_queue.scheduleEvent(event_queue.timeNow(), Event::make_blossom({v, u, e}));
    } else {
        assert(u_blossom->label == even);

        reducedGraph.forEdgesOf(u, [this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            if (edge_in_matching[id] || vertex_path[v] != path_root || 
                search_shell[v] != starting_shell ||
                union_find.find(u) == union_find.find(v)) return;
            
            // #if DEBUG_LOGGING
            // std::cerr << "CONSIDER NEIGHBOR " << v << std::endl;
            // #endif

            auto edge_slack = slack({u, v, id});
            auto v_blossom = split_find_min.list(v)->id;

            if (v_blossom->label == even) {
                // assert(edge_slack % 2 == 0);

                event_queue.scheduleEvent(
                    event_queue.timeNow() + edge_slack / 2, 
                    Event::make_blossom({u, v, id})
                );
            } else if (v_blossom->label == odd) { 
                auto [key, edge] = split_find_min.currentKey(v);
                if (edge == no_edge || edge_slack < slack(edge)) {
                    split_find_min.decreaseKey(
                        v, edge_slack + (v_blossom->t_odd - v_blossom->Delta), {u, v, id});
                }
            } else if (v_blossom->label == free) {
                auto min_slack = split_find_min.findMin(v_blossom->list).first - 
                                 (event_queue.timeNow() - v_blossom->Delta);
                
                split_find_min.decreaseKey(
                        v, edge_slack + (event_queue.timeNow() - v_blossom->Delta), {u, v, id});

                // #if DEBUG_LOGGING
                // auto min_edge =  split_find_min.findMin(v_blossom->list).second;
                // std::cerr << "min_slack " << min_slack << std::endl;
                // std::cerr << "min_edge " << min_edge << std::endl;
                // #endif
                
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

    for (auto it = B->subblossoms.begin(); std::next(it) != B->subblossoms.end(); it ++) {
        auto this_blossom = it->first;
        auto next_blossom = std::next(it)->first;
        auto [L1, L2] = split_find_min.split(this_blossom->base, this_blossom, next_blossom);
        this_blossom->list = L1;
        next_blossom->list = L2;

        // #if DEBUG_LOGGING
        // std::cerr << "LIST(" << this_blossom << ") = " << L1 << std::endl;
        // std::cerr << "LIST(" << next_blossom << ") = " << L2 << std::endl;
        // #endif
    }

    std::vector<NetworKit::node> to_schedule;

    for (auto [b, e] : B->subblossoms) {
        add_blossom(b, active_shell);
        b->parent = nullptr;
        b->t_root = event_queue.timeNow();
        b->Delta = B->Delta + (event_queue.timeNow() - B->t_odd);
        if (b->label == odd ) {
            b->t_odd  = event_queue.timeNow();
            if (!b->is_trivial())
                event_queue.scheduleEvent(event_queue.timeNow() + b->z0 / 2, Event::make_dissolve(b));
        } else if (b->label == even) {
            b->t_even = event_queue.timeNow();
            b->for_nodes([this, b, &to_schedule] (NetworKit::node v) {
                if (v != b->base)
                    union_find.link(b->base, v, b);
                to_schedule.push_back(v);
            });
        } else {
            auto [key, edge] = split_find_min.findMin(b->list);
            if (edge != no_edge) {
                auto [u, v, id] = edge;
                auto slack = key - (event_queue.timeNow() - b->Delta);
                event_queue.scheduleEvent(event_queue.timeNow() + slack, Event::make_grow(v, edge));
            }
        }
    }

    #if DEBUG_LOGGING
    print_search_state();
    #endif

    for (auto v : to_schedule) {
        schedule(v);
    }

    delete_blossom(B, active_shell);
    delete B;
}

std::list<GabowScalingMatching::Blossom*> 
GabowScalingMatching::blossoms_containing(NetworKit::node u, Blossom* until) {
    // #if DEBUG_LOGGING
    // std::cerr << "find blossoms on path to " << u << " until blossom ";
    // until->nodes_print(); std::cerr << std::endl;
    // #endif

    std::list<Blossom*> res;
    Blossom* iter = trivial_blossom[u];
    while (iter != until) { 
        res.push_front(iter);
        iter = iter->parent;
    }
    return res;
}

std::pair<std::list<std::pair<GabowScalingMatching::Blossom*, GabowScalingMatching::EdgeInfo>>, 
          std::list<std::pair<GabowScalingMatching::Blossom*, GabowScalingMatching::EdgeInfo>>>
GabowScalingMatching::split_subblossoms(
        std::list<std::pair<Blossom*, EdgeInfo>> sub_blossoms,
        Blossom* blossom) {
            
    auto iter = sub_blossoms.begin();
    while (iter->first != blossom) iter ++;
    return {{sub_blossoms.begin(), std::next(iter)}, {std::next(iter), sub_blossoms.end()}};
}

void GabowScalingMatching::blossom(EdgeInfo edge) {
    auto [u, v, id] = edge;

    if (union_find.find(u) == union_find.find(v)) return;

    #if DEBUG_LOGGING
    assert(slack(edge) == 0);
    #endif

    Blossom* u_blossom = get_blossom(u);
    Blossom* v_blossom = get_blossom(v);

    backtrack(u_blossom, v_blossom, edge);
}

void GabowScalingMatching::dissolveShell(OldBlossom* S) {
    S->dissolved = true;
    S->z = 0;

    auto distributed = event_queue.timeNow() - active_shell->t_active;
    active_shell->z -= 2 * distributed;
    active_shell_initial_dual -= 2 * distributed;
    add_distribution(active_shell, distributed);
    active_shell->t_active = event_queue.timeNow();
    
    if (S == active_shell) {
        #if DEBUG_LOGGING
        std::cerr << "DISSOLVING ACTIVE SHELL\n";
        #endif

        auto outer_shell = S->heavy_path_parent;

        std::vector<EdgeInfo> edges_to_check;

        if (S == highest_undissolved) {
            #if DEBUG_LOGGING
            std::cerr << "FINISH - OUTERMOST SHELL DISSOLVED\n";
            #endif

            finish_shell_search();
            if (S != path_root) {
                assert(outer_shell == path_root);
            }
            while (highest_undissolved != nullptr && highest_undissolved->dissolved)
                highest_undissolved = highest_undissolved->heavy_child;
        } else {
            assert(outer_shell != nullptr);

            if (outer_shell->searched) {
                #if DEBUG_LOGGING
                std::cerr << "FINISH - OUTER SHELL ALREADY SEARCHED\n";
                #endif
                finish_shell_search();
            } else {
                for (auto v : outer_shell->nodes) {
                    t_shell[v] = event_queue.timeNow();
                    Delta[v] = distribution_so_far(S->shell_index);
                    search_shell[v] = starting_shell;
                }
                for (auto b : outer_shell->shell_blossoms) {
                    b->list = split_find_min.init(b->node_list(), b);
                    b->Delta = 0;
                    b->label = free;
                }
                for (auto v : outer_shell->nodes) {
                    if (matched_vertex[v] == NetworKit::none)
                        event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, no_edge));
                    reducedGraph.forEdgesOf(v, 
                        [this, &edges_to_check] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
                            if (vertex_path[u] != path_root ||
                                search_shell[u] != starting_shell) return;
                            if (split_find_min.list(u)->id->label == even) {
                                edges_to_check.push_back({u, v, id});
                            }
                        }
                    );
                }
            }
        }
        if (outer_shell == nullptr) return;

        outer_shell->nodes.splice(outer_shell->nodes.end(), S->nodes);
        outer_shell->shell_blossoms.splice(outer_shell->shell_blossoms.end(), std::move(S->shell_blossoms));  
        outer_shell->heavy_child = active_shell->heavy_child;
        if (outer_shell->heavy_child != nullptr)
            outer_shell->heavy_child->heavy_path_parent = outer_shell;
        active_shell = outer_shell;
        active_shell_initial_dual = current_shell_duals[active_shell->shell_index] - 
                                    2 * distribution_so_far(active_shell->shell_index + 1);
        active_shell->t_active = event_queue.timeNow();
        active_shell->searched = true;
        if (search_done) return;
        if (active_shell != old_root) {
            event_queue.scheduleEvent(
                event_queue.timeNow() + active_shell->z / 2,
                Event::make_dissolveShell(active_shell)
            );
        } else {
            t_undissolvable = event_queue.timeNow();
        }

        for (auto e : edges_to_check)
            event_queue.scheduleEvent(
                event_queue.timeNow() + slack(e), 
                Event::make_grow(e.v, e)
            );
    } else {
        assert(S == active_shell->heavy_child);
        #if DEBUG_LOGGING
        std::cerr << "DISSOLVING INNER SHELL\n";
        #endif

        add_distribution(S, event_queue.timeNow() - S->t_inner);

        std::vector<EdgeInfo> edges_to_check;

        if (S->searched) {
            #if DEBUG_LOGGING
            std::cerr << "FINISH - INNER SHELL ALREADY SEARCHED\n";
            #endif
            S->t_inner = event_queue.timeNow();
            finish_shell_search();
        } else {
            for (auto v : S->nodes) {
                t_shell[v] = event_queue.timeNow();
                Delta[v] = distribution_so_far(S->shell_index + 1);
                search_shell[v] = starting_shell;
            };
            for (auto b : S->shell_blossoms) {
                b->list = split_find_min.init(b->node_list(), b);
                b->Delta = 0;
                b->label = free;
            }
            for (auto v : S->nodes) {
                if (matched_vertex[v] == NetworKit::none)
                    event_queue.scheduleEvent(event_queue.timeNow(), Event::make_grow(v, no_edge));
                reducedGraph.forEdgesOf(v, 
                    [this, &edges_to_check] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
                        if (vertex_path[u] != path_root ||
                            search_shell[u] != starting_shell) return;
                        if (split_find_min.list(u)->id->label == even) {
                            edges_to_check.push_back({u, v, id});
                        }
                    }
                );
            }
        }

        active_shell->nodes.splice(active_shell->nodes.end(), S->nodes);
        active_shell->shell_blossoms.splice(active_shell->shell_blossoms.end(), std::move(S->shell_blossoms));
        active_shell->heavy_child = S->heavy_child;
        if (active_shell->heavy_child != nullptr) {
            active_shell->heavy_child->heavy_path_parent = active_shell;
            active_shell->heavy_child->t_inner = event_queue.timeNow();
            if (search_done) return;
            event_queue.scheduleEvent(
                event_queue.timeNow() + active_shell->heavy_child->z / 2,
                Event::make_dissolveShell(active_shell->heavy_child)
            );
        }
        if (search_done) return;
        for (auto e : edges_to_check)
            event_queue.scheduleEvent(
                event_queue.timeNow() + slack(e), 
                Event::make_grow(e.v, e)
            );
    }

    #if DEBUG_LOGGING
    active_shell->nodes.sort();
    #endif
}

bool GabowScalingMatching::backtrack(Blossom* u, Blossom* v, EdgeInfo edge) {
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

            for (auto [blossom, e] : u_path) blossom->visited = false; 
            for (auto [blossom, e] : v_path) blossom->visited = false; 
            u->visited = false; v->visited = false;

            create_new_blossom(u, v, edge, u_path, v_path);

            return false;
        }
    }

    for (auto [blossom, e] : u_path) blossom->visited = false;
    for (auto [blossom, e] : v_path) blossom->visited = false;
    u->visited = false; v->visited = false;

    // The two paths don't intersect - an augmenting path was found ending the stage
    #if DEBUG_LOGGING
    std::cerr << "FINISH - AUGMENTING PATH FOUND\n";
    print_backtrack(u, v, edge, u_path, v_path);
    #endif
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
        Blossom* u, Blossom* v, EdgeInfo edge, 
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    
    cut_path_at(u_path, v_path.size() > 0 ? v_path.back().blossom : v, v);
    cut_path_at(v_path, u_path.size() > 0 ? u_path.back().blossom : u, u);
    if (v_path.size() > 0 && v_path.back().blossom == u) u_path.clear();
    if (u_path.size() > 0 && u_path.back().blossom == v) v_path.clear();

    #if DEBUG_LOGGING
    std::cerr << "Creating new blossom" << std::endl;
    print_backtrack(u, v, edge, u_path, v_path);
    #endif

    Blossom* base = v_path.size() > 0 ? v_path.back().blossom : v;
    std::list<std::pair<Blossom*, EdgeInfo>> subblossoms;
    
    for (int i = v_path.size() - 1; i > 0; -- i) {
        subblossoms.emplace_back(v_path[i-1].blossom, v_path[i].edge);
    }

    if (v_path.size() > 0) 
        subblossoms.emplace_back(v, v_path[0].edge);
    subblossoms.emplace_back(u, reverse(edge));

    for (int i = 0; i < u_path.size(); ++ i) {
        subblossoms.emplace_back(u_path[i].blossom, reverse(u_path[i].edge));
    }

    #if DEBUG_LOGGING
    std::cerr << "SUBBLOSSOMS:\n";
    for (auto [B, e] : subblossoms) {
        B->short_print();
        std::cerr << " " << e << std::endl;
    }
    #endif

    Blossom* new_blossom = new Blossom {
        base->base, base->base, // base, initial_base
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
            b->t_even = event_queue.timeNow();
            b->Delta = b->Delta + (event_queue.timeNow() - b->t_odd);
            b->for_nodes([this, b, new_blossom, &to_schedule] (NetworKit::node v) {
                to_schedule.push_back(v);
                if (v != b->base)
                    union_find.link(b->base, v, new_blossom);
            });
        }
        union_find.link(new_blossom->base, b->base, new_blossom);
        delete_blossom(b, active_shell);
    }
    add_blossom(new_blossom, active_shell);

    for (auto v : to_schedule) schedule(v);

    #if DEBUG_LOGGING
    std::cerr << "Created new blossom: " << std::endl;
    new_blossom->short_print(); std::cerr << std::endl;
    #endif
}

void GabowScalingMatching::cut_path_at(
        std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut2) {
    int index = -1;
    for (unsigned i = 0; i < path.size(); ++ i) {
        if (path[i].blossom == cut || path[i].blossom == cut2) {
            index = i;
            break;
        }
    }
    if (index != -1) path.resize(index + 1);
}

MaximumMatching::intedgeweight GabowScalingMatching::y(NetworKit::node v) {
    MaximumMatching::intedgeweight basic = y0[v] + Delta[v] + 
            std::max(0, std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    auto B = split_find_min.list(v)->id;

    if (B->label == odd) {
        return basic + B->Delta + event_queue.timeNow() - B->t_odd;
    } else if (B->label == even) {
        return basic + B->Delta - (event_queue.timeNow() - B->t_even);
    }
    return basic + B->Delta;
}

MaximumMatching::intedgeweight GabowScalingMatching::z(Blossom* B) {
    if (B->label == odd) {
        return B->z0 - 2 * (event_queue.timeNow() - B->t_odd);
    } else if (B->label == even) {
        return B->z0 + 2 * (event_queue.timeNow() - B->t_even);
    } else {
        return B->z0;
    }
}

MaximumMatching::intedgeweight GabowScalingMatching::shell_z() {
    return active_shell == old_root ? 0 :
        (active_shell_initial_dual - 2 * (event_queue.timeNow() - active_shell->t_active));
}

MaximumMatching::intedgeweight GabowScalingMatching::slack(EdgeInfo e) {
    auto [u, v, id] = e;
    auto u_B = get_blossom(u);
    auto v_B = get_blossom(v); 
    return y(u) + y(v) - current_w[id] + (u_B == v_B ? z(u_B) : 0) + outer_blossom_dual + shell_z();
}

bool GabowScalingMatching::is_exposed(Blossom *b) {
    return matched_vertex[b->base] == NetworKit::none;
}

bool GabowScalingMatching::is_even(NetworKit::node v) {
    return split_find_min.list(v)->id->label == even;
}

bool GabowScalingMatching::is_odd(NetworKit::node v) {
    return split_find_min.list(v)->id->label == odd;
}

void GabowScalingMatching::add_distribution(OldBlossom* S, MaximumMatching::intedgeweight distribution) {
    shell_distribution.add(S->shell_index, distribution);
}

MaximumMatching::intedgeweight GabowScalingMatching::distribution_so_far(int shell_index) {
    return shell_distribution.sum(shell_index - 1);
}

GabowScalingMatching::Blossom* GabowScalingMatching::get_blossom(NetworKit::node v) {
    auto B = split_find_min.list(v)->id;
    if (B->label == even) {
        return union_find.find(B->base);
    } else {
        return B;
    }
}

bool GabowScalingMatching::OldBlossom::is_heavy_path_root() {
    return heavy_path_parent == nullptr;
}

void GabowScalingMatching::OldBlossom::for_nodes(const std::function<void(NetworKit::node)>& handle) {
    for (auto n : nodes)
        handle(n);
    for (auto b : subblossoms)
        b->for_nodes(handle);
    if (heavy_child != nullptr) 
        heavy_child->for_nodes(handle);
}

void GabowScalingMatching::OldBlossom::for_blossoms(const std::function<void(OldBlossom*)>& handle) {
    handle(this);
    for (auto b : subblossoms)
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

GabowScalingMatching::EdgeInfo GabowScalingMatching::reverse(const EdgeInfo& info) {
    auto [u, v, id] = info;
    return { v, u, id };
}

bool GabowScalingMatching::EdgeInfo::operator==(const GabowScalingMatching::EdgeInfo& other) const {
    return u == other.u && v == other.v && id == other.id;
}

bool GabowScalingMatching::EdgeInfo::operator!=(const GabowScalingMatching::EdgeInfo& other) const {
    return !(*this == other);
}

bool GabowScalingMatching::EdgeInfo::operator<(const GabowScalingMatching::EdgeInfo& other) const {
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
    for (auto b : subblossoms) {
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

std::ostream& operator<<(std::ostream &out, const GabowScalingMatching::EdgeInfo& edge) {
    out << "(" << node_to_str(edge.u) << ", " << node_to_str(edge.v) << ")";
    return out;
}

std::ostream& operator<<(std::ostream &out, const GabowScalingMatching::Event& event) {
    switch (event.type) {
        case GabowScalingMatching::Event::Type::grow:
            out << "grow(" << event.args.v << ", " << event.args.e << ")";
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
        Blossom* u, Blossom* v, EdgeInfo edge, 
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
    for (int i = 0; i < shells.size(); ++ i) std::cerr << distribution_so_far(i + 1) << " ";
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

    MaximumMatching::intedgeweight basic = y0[v] + Delta[v] + 
            std::max(0, std::min(event_queue.timeNow(), t_undissolvable) - t_shell[v]);
    auto B = split_find_min.list(v)->id;
    
    MaximumMatching::intedgeweight y;
    
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

std::string GabowScalingMatching::label_to_str(BlossomLabel label) {
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