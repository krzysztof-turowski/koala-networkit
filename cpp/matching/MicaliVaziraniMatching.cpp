#include <matching/MaximumMatching.hpp>

namespace Koala {

MaximumCardinalityMatching::MaximumCardinalityMatching(NetworKit::Graph &graph): graph(graph) { }

const std::map<NetworKit::node, NetworKit::node>& MaximumCardinalityMatching::getMatching() const {
    assureFinished();
    return matching;
}

MicaliVaziraniMatching::MicaliVaziraniMatching(NetworKit::Graph &graph):
        MaximumCardinalityMatching(graph),
        V(graph.upperNodeIdBound()),
        E(graph.upperEdgeIdBound()),
        candidates(graph.upperNodeIdBound()),
        bridges(2 * graph.upperNodeIdBound() + 1),
        bloom_bases(graph.upperNodeIdBound()) {

        }

template<class C>
void print_nodes(const C& P) {
    for (auto v : P) std::cerr << v << " ";
    std::cerr << std::endl;
}

void print_nodes(const std::vector<std::pair<NetworKit::node, NetworKit::edgeid>>& P) {
    for (auto [v, e] : P) std::cerr << v << " ";
    std::cerr << std::endl;
}

std::string level_to_str(int level) {
    return level == 1e9 ? "inf" : std::to_string(level);
}

std::string node_to_str(NetworKit::node vertex) {
    return vertex == NetworKit::none ? "-" : std::to_string(vertex);
}

void MicaliVaziraniMatching::run() {
    #if DEBUG_LOGGING 
    std::cerr << "GRAPH:\n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << " : ";
        graph.forEdgesOf(v, [] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
            std::cerr << u << " "; // "(" << id << ") ";
        });
        std::cerr << std::endl;
    });
    std::cerr << std::endl;
    graph.forEdges([] (NetworKit::node u, NetworKit::node v) {
        std::cerr << u << " " << v << std::endl;
    });
    std::cerr << std::endl;
    #endif

    color_counter = 1;

    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        E[static_cast<unsigned int>(id)].u = u;
        E[static_cast<unsigned int>(id)].v = v;
        // E[static_cast<unsigned int>(id)].matched = false;
    });

    graph.forNodes([this] (NetworKit::node vertex) {
        V[vertex].match = NetworKit::none;
        V[vertex].match_edge = NetworKit::none;
    });

    do {
        augmentation_happened = false;
        color_counter = 1;

        graph.forNodes([this] (NetworKit::node vertex) {
            V[vertex].parent = NetworKit::none;
            V[vertex].parent_edge = NetworKit::none;
            V[vertex].even_level = inf_level;
            V[vertex].odd_level = inf_level;
            V[vertex].bloom = nullptr;
            V[vertex].predecessors.clear();
            V[vertex].pred_it = 0;
            V[vertex].successors.clear();
            V[vertex].children.clear();
            V[vertex].count = 0;
            V[vertex].erased = false;
            V[vertex].visited = false;
            V[vertex].color = no_color;
        });

        graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            id = static_cast<unsigned int>(id);
            E[id].visited = false;
            E[id].type = EdgeData::Type::none;
        });

        for (int i = 0; i < graph.upperNodeIdBound(); ++ i) {
            candidates[i].clear();
            bridges[2 * i + 1].clear();
        }

        bloom_bases.reset();
        erase_queue.clear();

        search();

        #if DEBUG_LOGGING
        check_consistency();
        #endif

        for (auto B : current_blooms) delete B;
        current_blooms.clear();
    } while (augmentation_happened);

    graph.forNodes([this] (NetworKit::node vertex) {
        matching[vertex] = V[vertex].match;
    });
    hasRun = true;
}

void MicaliVaziraniMatching::search() {
    #if DEBUG_LOGGING
    std::cerr << "=============================\n";
    std::cerr << "SEARCH" << std::endl;
    #endif

    graph.forNodes([this] (NetworKit::node vertex) {
        if (exposed(vertex)) {
            set_level(vertex, 0);
        }
    });

    iter = max_iter = 0;
    augmentation_happened = false;
    for (iter = 0; iter < candidates.size() && !augmentation_happened; iter ++) {
        #if DEBUG_LOGGING
        if (candidates[iter].size() > 0) {
            std::cerr << "------------------------------------------------------\n";
            std::cerr << "LEVEL " << iter << std::endl;
            print_state();
        }
        #endif

        for (auto v : candidates[iter]) {
            #if DEBUG_LOGGING
            std::cerr << "CHECK " << v << std::endl;
            #endif

            graph.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid e) {
                e = static_cast<int>(e);

                if (V[u].erased || E[e].type != EdgeData::Type::none ||
                    (iter % 2 == 0 && u == V[v].match) || (iter % 2 == 1 && u != V[v].match))
                    return;

                #if DEBUG_LOGGING
                std::cerr << "LOOK AT " << u << std::endl;
                #endif

                if (min_level(u) >= iter + 1) {
                    E[e].type = EdgeData::Type::prop;
                    if (min_level(u) > iter + 1) {
                        set_level(u, iter + 1);
                    }
                    V[u].predecessors.push_back(v);
                    V[u].count ++;
                    V[v].successors.push_back(u);

                    #if DEBUG_LOGGING
                    std::cerr << "MAKE " << v << " A PRED OF " <<  u << std::endl;
                    #endif
                } else {
                    E[e].type = EdgeData::Type::bridge;
                    auto t = tenacity(u, v);
                    if (t < inf_level) {
                        bridges[t].push_back(e);

                        #if DEBUG_LOGGING
                        std::cerr << "BRIDGE FOUND " << u << " " <<  v << " WITH TENACITY " << t << " AND ID " << e << std::endl;
                        #endif
                    }
                }
            });
        }

        #if DEBUG_LOGGING
        print_state();
        #endif 

        for (int i = 0; i < bridges[2 * iter + 1].size(); ++ i) {
            auto e = bridges[2 * iter + 1][i];
            auto s = E[e].v;
            auto t = E[e].u;
            if (!V[s].erased && !V[t].erased)
                bloss_aug(s, t, e);
        }
    }
}

void MicaliVaziraniMatching::flip(NetworKit::node u, NetworKit::node v) {
    #if DEBUG_LOGGING
    std::cerr << "FLIP " << u << " " << v << std::endl;
    #endif
    
    if(V[u].erased || V[v].erased || V[u].match == v) return;

    V[u].match = v;
    V[v].match = u;

    V[u].erased = true;
    V[v].erased = true;
    erase_queue.push_back(u);
    erase_queue.push_back(v);
}

bool MicaliVaziraniMatching::openingDfs(NetworKit::node cur, NetworKit::node bcur, NetworKit::node b) {
    #if DEBUG_LOGGING
    std::cerr << "OPENING DFS " << cur << " " << bcur << " " << b << std::endl;
    #endif
    
    if(bcur == b) {
        augumentPath(cur, bcur);
        return true;
    }

    #if DEBUG_LOGGING
    for(auto a: V[bcur].children)
        std::cerr << a.first << " " << a.second << " " << V[a.second].color << " " << V[bcur].color << std::endl;
    #endif

    for(auto a: V[bcur].children) { 
        if((a.second == b || V[a.second].color == V[bcur].color) && openingDfs(a.first,a.second,b)) {
            augumentPath(cur, bcur);
            flip(bcur,a.first);
            return true;
        }
    }
    return false;
}

void MicaliVaziraniMatching::augumentPath(NetworKit::node u, NetworKit::node v, bool initial) {
    #if DEBUG_LOGGING
    std::cerr << "AUGMENT PATH " << u << " " << v << " " << initial << std::endl;
    #endif
    
    if(u == v) return;
    if(!initial && outer(u)) {
        int x = V[u].predecessors.front();

        int idx = 0;
        while (base_star(V[x].predecessors[idx]) != base_star(x)) {
            idx ++; 
        }

        u = V[x].predecessors[idx];
        flip(x,u);

        augumentPath(u,v);
    } else {
        auto [u3, v3, u2, v2] = get_bridge(u);

        flip(u3, v3);
        bool openingDfsSucceed1 = openingDfs(u3,u2,u);

        int v4 = base(u);
        bool openingDfsSucceed2 = openingDfs(v3,v2,v4);
        augumentPath(v4, v);
    }
}

void MicaliVaziraniMatching::bloss_aug(NetworKit::node g, NetworKit::node r, NetworKit::edgeid e) {    
    if (V[g].bloom == V[r].bloom && V[g].bloom != nullptr) return;

    #if DEBUG_LOGGING
    std::cerr << "\nBLOSS AUG " << g << " " <<  r << std::endl;
    #endif

    auto v_G = (V[g].bloom == nullptr) ? g : base_star(V[g].bloom);
    auto v_R = (V[r].bloom == nullptr) ? r : base_star(V[r].bloom);
    auto green_root = v_G;
    auto red_root = v_R;

    if (v_G == v_R) return;

    int green_color = ++ color_counter;
    int red_color = ++ color_counter;
    V[v_G].color = green_color;
    V[v_G].parent = v_G;
    V[v_R].color = red_color;
    V[v_R].parent = r;
    auto barrier = v_G;
    bridge_support = { v_G, v_R };
    bloom_found = false;

    #if DEBUG_LOGGING
    std::cerr << "MARK " << v_G << " green " << green_color << std::endl;
    std::cerr << "MARK " << v_R << " red " << red_color << std::endl;
    #endif

    while ((!exposed(v_G) || !exposed(v_R)) && !bloom_found) {
        #if DEBUG_LOGGING
        std::cerr << "DFS " << node_to_str(v_R) << " " <<  node_to_str(v_G) << " " << node_to_str(barrier) << std::endl;
        #endif

        // if (v_G == NetworKit::none || v_R == NetworKit::none) return;

        if (min_level(v_R) >= min_level(v_G)) {
            #if DEBUG_LOGGING
            std::cerr << "ADVANCE RED\n";
            #endif

            red_dfs_step(v_R, red_color, v_G, green_color, r, barrier);
        } else {
            #if DEBUG_LOGGING
            std::cerr << "ADVANCE GREEN\n";
            #endif

            green_dfs_step(v_G, green_color, v_R, red_color, r, barrier);
        }
    }

    Bloom *B = new Bloom { v_G, green_color, red_color, g, green_root, r, red_root };
    current_blooms.push_back(B);
    V[B->base].color = no_color;

    #if DEBUG_LOGGING
    std::cerr << "BLOOM FOUND " << std::endl;
    std::cerr << "ADDRESS = " << B << std::endl;
    std::cerr << "BASE = " << B->base << std::endl;
    std::cerr << "GREEN COLOR = " << B->green_color << std::endl;
    std::cerr << "RED COLOR = " << B->red_color << std::endl;
    std::cerr << "GREEN PEAK = " << B->green_peak << std::endl;
    std::cerr << "RED PEAK = " << B->red_peak << std::endl;
    std::cerr << "NODES: ";
    for (auto y : bridge_support)
        if (V[y].color != no_color) 
            std::cerr << y << " ";
    std::cerr << std::endl;
    #endif

    for (auto y : bridge_support) {
        if (V[y].color == no_color) continue;

        V[y].bloom = B;
        bloom_bases.link(B->base, y);
        set_level(y, 2 * iter + 1 - min_level(y));
        
        if (inner(y)) {
            graph.forEdgesOf(y, [this] (NetworKit::node y, NetworKit::node u, NetworKit::edgeid e) {
                e = static_cast<int>(e);
                if (E[e].type != EdgeData::Type::bridge) return;

                int t = tenacity(y, u);
                if (t < inf_level) 
                    bridges[t].push_back(e);
            });
        }
    }

    if (bloom_found) return;

    augmentation_happened = true;
    augumentPath(v_R, v_G, true);

    erase(erase_queue);
}

void MicaliVaziraniMatching::red_dfs_step(
        NetworKit::node& v_R, int red_color, NetworKit::node& v_G, int green_color, 
        NetworKit::node r, NetworKit::node& barrier) {
    
    while (V[v_R].pred_it < V[v_R].predecessors.size()) {
        auto u = V[v_R].predecessors[V[v_R].pred_it ++];
        auto x = u;

        if (V[u].erased) continue;

        if (V[u].bloom != nullptr) u = base_star(V[u].bloom);

        if (V[u].color == no_color) {
            V[v_R].children.push_back({x, u});
            V[u].color = red_color;
            V[u].parent = v_R;

            v_R = u;
            bridge_support.push_back(v_R);

            return;
        } else if (u == v_G) {
            V[v_R].children.push_back({x, u});

            if (u != barrier) {
                // Green dfs is not at the barrier
                // Take over it's center and force it to backtrack

                v_G = V[v_G].parent;

                V[u].color = red_color;
                V[u].parent = v_R;
                v_R = u;

                return;
            }
        }
    }

    if (v_R ==  r) {
        // Backtracked to root - a bottleneck was found
        bloom_found = true;
    } else {
        // Backtrack to parent
        v_R = V[v_R].parent;
    }
}

void MicaliVaziraniMatching::green_dfs_step(
        NetworKit::node& v_G, int green_color, NetworKit::node& v_R, int red_color, 
        NetworKit::node r, NetworKit::node& barrier) {
    
    while (V[v_G].pred_it < V[v_G].predecessors.size()) {
        auto u = V[v_G].predecessors[V[v_G].pred_it ++];
        auto x = u;

        if (V[u].erased) continue;

        if (V[u].bloom != nullptr) u = base_star(V[u].bloom);

        if (V[u].color == no_color) {
            V[v_G].children.push_back({x, u});
            V[u].color = green_color;
            V[u].parent = v_G;

            v_G = u;
            bridge_support.push_back(v_G);

            return;
        } else if (u == v_R) {
            V[v_G].children.push_back({x, u});
        }
    }

    if (v_G == barrier) {
        // Backtracked to the barrier
        // Take over the center of red dfs and force it to backtrack

        barrier = v_R;
        v_G = v_R;
        V[v_G].color = green_color;

        v_R = V[v_R].parent;
    } else {
        // Backtrack
        v_G = V[v_G].parent;
    }
}

// void MicaliVaziraniMatching::dfs_step(
//         NetworKit::node& v_1, int color_1, NetworKit::node& v_2, int color_2, 
//         int red_color, int r, NetworKit::node& barrier) {
    
//     while (V[v_1].pred_it < V[v_1].predecessors.size()) {
//         auto u = V[v_1].predecessors[V[v_1].pred_it ++];
//         auto x = u;

//         #if DEBUG_LOGGING
//         std::cerr << "LOOK AT " << u << " PRED OF " << v_1 << std::endl;;
//         #endif

//         if (V[u].erased) continue;

//         if (V[u].bloom != nullptr) u = base_star(V[u].bloom);

//         if (V[u].color == no_color) {
//             V[v_1].children.push_back({x, u});
//             V[u].color = color_1;
//             V[u].parent = v_1;

//             #if DEBUG_LOGGING
//             std::cerr << "PARENT(" << u << ") = " << v_1 << std::endl;
//             std::cerr << "MARK " << u << " WITH " << color_1 << std::endl;
//             #endif

//             v_1 = u;
//             bridge_support.push_back(v_1);

//             #if DEBUG_LOGGING
//             std::cerr << "MOVE TO " << v_1 << "\n";
//             #endif

//             return;
//         } else if (u == v_2) {
//             V[v_1].children.push_back({x, u});
//             if (color_1 == red_color && u != barrier) {
//                 v_2 = V[v_2].parent;

//                 #if DEBUG_LOGGING
//                 std::cerr << "BACKTRACK GREEN TO " << v_2 << "\n";
//                 #endif

//                 V[u].color = color_1;
//                 V[u].parent = v_1;

//                 #if DEBUG_LOGGING
//                 std::cerr << "PARENT(" << u << ") = " << v_1 << std::endl;
//                 std::cerr << "MARK " << u << " WITH " << color_1 << std::endl;
//                 #endif

//                 v_1 = u;

//                 #if DEBUG_LOGGING
//                 std::cerr << "TAKE OVER " << v_1 << "\n";
//                 #endif
//                 return;
//             }
//         }
//     }

//     if (color_1 == red_color) {
//         if (v_1 ==  r) {
//             bloom_found = true;

//             #if DEBUG_LOGGING
//             std::cerr << "BLOOM FOUND\n";
//             #endif
//         } else {
//             v_1 = V[v_1].parent;

//             #if DEBUG_LOGGING
//             std::cerr << "BACKTRACK TO " << v_1 << "\n";
//             #endif
//         }
//     } else {
//         if (v_1 == barrier) {
//             barrier = v_2;
//             v_1 = v_2;
//             V[v_1].color = color_1;
//             v_2 = V[v_2].parent;

//             #if DEBUG_LOGGING
//             std::cerr << "GIVE OVER " << v_1 << " TO GREEN\n";
//             std::cerr << "MARK " << v_1 << " WITH " << color_1 << std::endl;
//             #endif
//         } else {
//             v_1 = V[v_1].parent;

//             #if DEBUG_LOGGING
//             std::cerr << "BACKTRACK TO " << v_1 << "\n";
//             #endif
//         }
//     }
// }

void MicaliVaziraniMatching::erase(std::vector<NetworKit::node>& Y) {
    for (int i = 0; i < Y.size(); ++ i) {
        auto y = Y[i];

        for (auto z : V[y].successors) {
            if (V[z].erased) continue;

            V[z].count --;
            if (V[z].count == 0) {
                V[z].erased = true;
                Y.push_back(z);
            }
        }
    }
}

// std::list<NetworKit::node> MicaliVaziraniMatching::find_path(
//     NetworKit::node high, NetworKit::node low, Bloom* B, int color) {

//     #if DEBUG_LOGGING
//     std::cerr << "FIND PATH " << high << " " << low << " in " << B << " MARKED " << color << std::endl;
//     #endif 

//     if (high == low) {
//         return {{ high }, {}};
//     }

//     auto v = high;
//     auto u = high;
//     NetworKit::edgeid ue;
//     while (u != low) {
//         #if DEBUG_LOGGING
//         std::cerr << std::endl;
//         std::cerr << "V = " << v << std::endl;
//         std::cerr << "U = " << u << std::endl;
//         #endif 

//         bool unvisited = false;
//         for (auto [q, e] : V[v].predecessors) {
//             if (E[static_cast<int>(e)].visited) continue;

//             #if DEBUG_LOGGING
//             std::cerr << "FOUND UNVISITED PRED " << q <<  "\n";
//             #endif

//             if (V[v].bloom == nullptr || V[v].bloom == B) {
//                 E[static_cast<int>(e)].visited = true;

//                 #if DEBUG_LOGGING
//                 std::cerr << "MARK (" << E[static_cast<int>(e)].u << ", " << E[static_cast<int>(e)].v <<  ") visited\n";
//                 #endif

//                 u = q;
//             } else {
//                 u = V[v].bloom->base;
//             }
            
//             if (u == low) break;
//         }

//         if (!unvisited && u != low) {
//             #if DEBUG_LOGGING
//             std::cerr << "BACK TO PARENT OF " << v << " = " << V[v].parent << std::endl;
//             #endif 

//             // V[v].visited = false;
//             v = V[v].parent;
//         }  
//     }

//     #if DEBUG_LOGGING
//     std::cerr << "PARENT(" << u << ") = " << v << std::endl;
//     #endif

//     std::list<NetworKit::node> P;
//     std::list<NetworKit::edgeid> EP;
//     V[u].parent = v;
//     V[u].parent_edge = ue;
//     while (u != high) {
//         #if DEBUG_LOGGING
//         std::cerr << "APPEND " << u << " WITH PARENT " << V[u].parent << std::endl;
//         #endif 

//         P.push_back(u);
//         EP.push_back(V[u].parent_edge);
//         u = V[u].parent;
//     } 
//     P.push_back(u);
//     std::reverse(P.begin(), P.end());
//     std::reverse(EP.begin(), EP.end());

//     #if DEBUG_LOGGING
//     std::cerr << "PREPATH:\n";
//     print_path(P, EP);
//     #endif 

//     auto it = P.begin();
//     auto eit = EP.begin();
//     while (std::next(it) != P.end() && it != P.end()) {
//         auto x = it;
//         auto y = std::next(it);
//         auto e = eit;
//         if (V[*x].bloom != nullptr && V[*x].bloom != B) {
//             #if DEBUG_LOGGING
//             std::cerr << "OPENING " << *x << " IN " << V[*x].bloom << " NOT IN " << B << std::endl;
//             #endif 

//             auto [O, EO] = open(*x);

//             #if DEBUG_LOGGING
//             std::cerr << "OPEN RESULT:\n";
//             print_path(O, EO);
//             #endif

//             P.splice(x, O);
//             EP.splice(e, EO);

//             it = std::prev(x);
//             eit = std::next(e);

//             P.erase(x);
//             P.erase(y);
//             EP.erase(e);

//             #if DEBUG_LOGGING
//             std::cerr << "PATH AFTER OPENING:\n";
//             print_path(P, EP);
//             std::cerr << "LOOKING AT " << *it << " AFTER " << *std::prev(it) << std::endl;
//             #endif 
//         } else {
//             it = std::next(it);
//             eit = std::next(eit);
//         }
//     }

//     return {P, EP};
// }

// std::list<NetworKit::node> MicaliVaziraniMatching::open(NetworKit::node x) {
//     #if DEBUG_LOGGING
//     std::cerr << "OPEN " << x << std::endl;
//     #endif

//     auto B = V[x].bloom;
//     auto b = B->base;
//     if (outer(x)) {
//         #if DEBUG_LOGGING
//         std::cerr << "OUTER" << std::endl;
//         #endif

//         return find_path(x, b, B, V[x].color);
//     } else {
//         auto gp = B->green_color;
//         auto rp = B->red_peak;
//         if (V[x].color == B->green_color) {
//             #if DEBUG_LOGGING
//             std::cerr << "MARKED LEFT" << std::endl;
//             #endif      

//             auto P = find_path(gp, x, B, B->green_color);
//             auto P2 = find_path(rp, b, B, B->red_color);
//             std::reverse(P.begin(), P.end());
//             P.splice(P.end(), P2);
//             return P;
//         } else if (V[x].color == B->red_color) {
//             #if DEBUG_LOGGING
//             std::cerr << "MARKED RIGHT" << std::endl;
//             #endif      

//             auto P = find_path(rp, x, B, B->red_color);
//             auto P2 = find_path(gp, b, B, B->green_color);
//             std::reverse(P.begin(), P.end());
//             P.splice(P.end(), P2);
//             return P;
//         }
//     }
// }

NetworKit::node MicaliVaziraniMatching::base_star(Bloom* bloom) {
    return bloom_bases.find(bloom->base);
}

NetworKit::node MicaliVaziraniMatching::base_star(NetworKit::node vertex) {
    return V[vertex].bloom == nullptr ? vertex : base_star(V[vertex].bloom);
}

NetworKit::node MicaliVaziraniMatching::base(NetworKit::node vertex) {
    return V[vertex].bloom == nullptr ? vertex : V[vertex].bloom->base;
}

bool MicaliVaziraniMatching::exposed(NetworKit::node vertex) {
    return V[vertex].match == NetworKit::none;
}

void MicaliVaziraniMatching::set_level(NetworKit::node vertex, int level) {
    if (level % 2 == 0) {
        V[vertex].even_level = level;
    } else {
        V[vertex].odd_level = level;
    }
    candidates[level].push_back(vertex);
}

int MicaliVaziraniMatching::min_level(NetworKit::node vertex) {
    return std::min(V[vertex].odd_level, V[vertex].even_level);
}

int MicaliVaziraniMatching::max_level(NetworKit::node vertex) {
    return std::max(V[vertex].odd_level, V[vertex].even_level);
}

int MicaliVaziraniMatching::tenacity(NetworKit::node u, NetworKit::node v) {
    return V[u].match == v ?
            V[u].odd_level  + V[v].odd_level  + 1 :
            V[u].even_level + V[v].even_level + 1;
}

bool MicaliVaziraniMatching::outer(NetworKit::node vertex) {
    return min_level(vertex) % 2 == 0;
}

bool MicaliVaziraniMatching::inner(NetworKit::node vertex) {
    return min_level(vertex) % 2 == 1;
}

std::tuple<NetworKit::node, NetworKit::node, NetworKit::node, NetworKit::node> 
MicaliVaziraniMatching::get_bridge(NetworKit::node vertex) {
    auto B = V[vertex].bloom;
    auto u3 = B->red_peak, v3 = B->green_peak, u2 = B->red_root, v2 = B->green_root;
    return V[vertex].color == B->red_color ?
        std::make_tuple(B->red_peak, B->green_peak, B->red_root, B->green_root) :
        std::make_tuple(B->green_peak, B->red_peak, B->green_root, B->red_root);
}

void MicaliVaziraniMatching::print_state() {
    graph.forNodes([this] (NetworKit::node vertex) {
        std::cerr << "VERTEX " << vertex << ": " << std::endl;
        std::cerr << "   MATCH       " << node_to_str(V[vertex].match) << std::endl;
        // std::cerr << "   MATCH EDGE  " << node_to_str(V[vertex].match_edge) << std::endl;
        std::cerr << "   EVEN LEVEL  " << level_to_str(V[vertex].even_level) << std::endl;
        std::cerr << "   ODD LEVEL   " << level_to_str(V[vertex].odd_level) << std::endl;
        // std::cerr << "   PARENT      " << node_to_str(V[vertex].parent) << std::endl;
        // std::cerr << "   PARENT EDGE " << node_to_str(V[vertex].parent_edge) << std::endl;
        std::cerr << "   PRED        "; print_nodes(V[vertex].predecessors);
        // std::cerr << "   SUCC        "; print_nodes(V[vertex].successors);
        std::cerr << "   BLOOM       " << V[vertex].bloom << std::endl;
        std::cerr << "   COLOR       " << V[vertex].color << std::endl;
        // std::cerr << "   EDGES       ";
        // graph.forEdgesOf(vertex, [this] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
            // std::cerr << u << "(" << id << ") ";
        // });
        // std::cerr << std::endl;
    });
    // for (auto B : current_blooms) {
    //     std::cerr << "BLOSSOM " << B << std::endl;
    //     std::cerr << "   BASE        " << B->base << std::endl;
    //     std::cerr << "   LEFT PEAK   " << B->left_peak << std::endl;
    //     std::cerr << "   RIGHT PEAK  " << B->right_peak << std::endl;
    //     std::cerr << "   BASE STAR   " << base_star(B) << std::endl;
    // }
}

void MicaliVaziraniMatching::check_consistency() {
    // std::cerr << "GRAPH:\n";
    // graph.forNodes([this] (NetworKit::node v) {
    //     std::cerr << v << " : ";
    //     graph.forEdgesOf(v, [] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid id) {
    //         std::cerr << u << "(" << id << ") ";
    //     });
    //     std::cerr << std::endl;
    // });
    // std::cerr << std::endl;
    // graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        
    // });
}

} /* namespace Koala */
