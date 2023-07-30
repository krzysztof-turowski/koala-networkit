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
        bridges(graph.upperNodeIdBound()) {}

template<class C>
void print_nodes(const C& P) {
    for (auto v : P) std::cerr << v << " ";
    std::cerr << std::endl;
}

void print_nodes(const std::vector<std::pair<NetworKit::node, NetworKit::edgeid>>& P) {
    for (auto [v, e] : P) std::cerr << v << " ";
    std::cerr << std::endl;
}

void MicaliVaziraniMatching::run() {
    #if DEBUG_LOGGING 
    std::cerr << "GRAPH:\n";
    graph.forNodes([this] (NetworKit::node v) {
        std::cerr << v << " : ";
        graph.forNeighborsOf(v, [] (NetworKit::node u) {
            std::cerr << u << " ";
        });
        std::cerr << std::endl;
    });
    std::cerr << std::endl;
    #endif

    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
        E[id].u = u;
        E[id].v = v;
        E[id].matched = false;
    });

    graph.forNodes([this] (NetworKit::node vertex) {
        V[vertex].match = NetworKit::none;
        V[vertex].match_edge = NetworKit::none;
    });

    do {
        augmentation_happened = false;

        graph.forNodes([this] (NetworKit::node vertex) {
            V[vertex].parent = NetworKit::none;
            V[vertex].parent_edge = NetworKit::none;
            V[vertex].even_level = std::numeric_limits<int>::max();
            V[vertex].odd_level = std::numeric_limits<int>::max();
            V[vertex].bloom = nullptr;
            V[vertex].predecessors.clear();
            V[vertex].successors.clear();
            V[vertex].anomalies.clear();
            V[vertex].count = 0;
            V[vertex].erased = false;
            V[vertex].visited = false;
            V[vertex].mark = VertexData::Mark::unmarked;
        });

        graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid id) {
            E[id].used = false;
            E[id].visited = false;
        });

        for (int i = 0; i < graph.upperNodeIdBound(); ++ i) {
            candidates[i].clear();
            bridges[i].clear();
        }

        bloom_bases.reset();

        search();

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
            V[vertex].even_level = 0;
            candidates[0].push_back(vertex);
        }
    });

    iter = 0;
    bool augmentation = false;
    while (!candidates[iter].empty() && !augmentation) {
        #if DEBUG_LOGGING
        std::cerr << "LEVEL " << iter << std::endl;
        #endif
        
        if (iter % 2 == 0) {
            for (auto v : candidates[iter]) {
                #if DEBUG_LOGGING
                std::cerr << "CHECK " << v << std::endl;
                #endif

                graph.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid e) {
                    if (V[u].erased || E[e].matched) return;

                    #if DEBUG_LOGGING
                    std::cerr << "LOOK AT NEIGHBOR " << u << std::endl;
                    #endif

                    if (V[u].even_level < std::numeric_limits<int>::max()) {
                        int j = (V[u].even_level + V[v].even_level) / 2;
                        bridges[j].push_back(e);

                        #if DEBUG_LOGGING
                        std::cerr << "BRIDGE FOUND " << u << " " <<  v << " AT LEVEL " << j << std::endl;
                        #endif
                    } else {
                        if (V[u].odd_level == std::numeric_limits<int>::max()) {
                            V[u].odd_level = iter + 1;
                            candidates[iter + 1].push_back(u);
                            
                            #if DEBUG_LOGGING
                            std::cerr << "QUEUE " << u << std::endl;
                            #endif
                        } 
                        if (V[u].odd_level == iter + 1) {
                            V[u].count ++;
                            V[u].predecessors.push_back({v, e});
                            V[v].successors.push_back(u);
                        } 
                        if (V[u].odd_level < iter) {
                            V[u].anomalies.push_back({v, e});
                        }
                    }
                });
            }
        }

        if (iter % 2 == 1) {
            for (auto v : candidates[iter]) {
                #if DEBUG_LOGGING
                std::cerr << "CHECK " << v << std::endl;
                #endif

                if (V[v].bloom != nullptr) continue;
                auto u = V[v].match;
                auto e = V[v].match_edge;

                #if DEBUG_LOGGING
                std::cerr << "LOOK AT MATCH " << u << std::endl;
                #endif

                if (V[u].odd_level < std::numeric_limits<int>::max()) {
                    int j = (V[u].odd_level + V[v].odd_level) / 2;
                    bridges[j].push_back(e);

                    #if DEBUG_LOGGING
                    std::cerr << "BRIDGE FOUND " << u << " " <<  v << " AT LEVEL " << j << std::endl;
                    #endif
                }
                if (V[u].even_level == std::numeric_limits<int>::max()) {
                    // V[u].predecessors.clear();
                    V[u].predecessors.push_back({v, V[u].match_edge});
                    // V[v].successors.clear();
                    V[v].successors.push_back(u);
                    V[u].count = 1;
                    V[u].even_level = iter + 1;
                    candidates[iter + 1].push_back(u);

                    #if DEBUG_LOGGING
                    std::cerr << "QUEUE " << u << std::endl;
                    #endif
                }
            }
        }

        #if DEBUG_LOGGING
        print_state();
        #endif 

        for (auto e : bridges[iter]) {
            auto s = E[e].u;
            auto t = E[e].v;
            if (!V[s].erased && !V[t].erased)
                bloss_aug(s, t, e);
        }

        iter ++;
    }
}


void MicaliVaziraniMatching::bloss_aug(NetworKit::node s, NetworKit::node t, NetworKit::edgeid e) {
    if (V[s].bloom == V[t].bloom && V[s].bloom != nullptr) return;

    #if DEBUG_LOGGING
    std::cerr << "BLOSS AUG " << s << " " <<  t << std::endl;
    #endif

    auto v_L = (V[s].bloom == nullptr) ? s : base_star(V[s].bloom);
    auto v_R = (V[t].bloom == nullptr) ? t : base_star(V[t].bloom);
    V[v_L].mark = VertexData::Mark::left;
    V[v_R].mark = VertexData::Mark::right;
    auto dcv = NetworKit::none;
    auto barrier = v_R;
    bloom_nodes = { v_L, v_R };
    bloom_found = false;

    while (!exposed(v_L) || !exposed(v_R)) {
        #if DEBUG_LOGGING
        std::cerr << "DFS " << v_L << " " <<  v_R << " " << dcv << " " << barrier << std::endl;
        #endif

        if (v_L == NetworKit::none || v_R == NetworKit::none) return;

        if (level(v_L) >= level(v_R)) {
            #if DEBUG_LOGGING
            std::cerr << "ADVANCE LEFT\n";
            #endif

            left_dfs(s, v_L, v_R, dcv, barrier);
        } else {
            #if DEBUG_LOGGING
            std::cerr << "ADVANCE RIGHT\n";
            #endif

            right_dfs(v_L, v_R, dcv, barrier);
        }
        if (bloom_found) {
            #if DEBUG_LOGGING
            std::cerr << "BLOOM FOUND " << std::endl;
            #endif
            if (dcv == NetworKit::none) return;

            Bloom *B = new Bloom { dcv, s, t, e };
            current_blooms.push_back(B);

            #if DEBUG_LOGGING
            std::cerr << "ADDRESS = " << B << std::endl;
            std::cerr << "BASE = " << B->base << std::endl;
            std::cerr << "LEFT PEAK = " << B->left_peak << std::endl;
            std::cerr << "RIGHT PEAK = " << B->right_peak << std::endl;
            std::cerr << "NODES: ";
            #endif

            V[dcv].mark = VertexData::Mark::unmarked;

            if (V[B->base].bloom == nullptr) {
                B->id = bloom_bases.create(B->base);
            } else {
                B->id = bloom_bases.link(V[B->base].bloom->id, B->base);
            }

            for (auto y : bloom_nodes) {
                if (V[y].mark == VertexData::Mark::unmarked) continue;

                #if DEBUG_LOGGING
                std::cerr << y << " ";
                #endif

                V[y].bloom = B;

                if (outer(y)) {
                    V[y].odd_level = 2 * iter + 1 - V[y].even_level;
                } else {
                    V[y].even_level = 2 * iter + 1  - V[y].odd_level;
                    candidates[V[y].even_level].push_back(y);
                    for (auto [z, e] : V[y].anomalies) {
                        int j = (V[y].even_level + V[z].even_level) / 2;
                        bridges[j].push_back(e);

                        #if DEBUG_LOGGING
                        std::cerr << "BRIDGE FOUND " << y << " " <<  z << " AT LEVEL " << j << std::endl;
                        #endif

                        E[e].used = true;
                    }
                }
            }

            #if DEBUG_LOGGING
            std::cerr << std::endl;
            #endif
            
            return;
        }
    }

    augmentation_happened = true;

    #if DEBUG_LOGGING
    std::cerr << "AUGMENT ON PATH " << v_L << " ~~ " << s << " - " << t << " ~~ " << v_R << std::endl;
    std::cerr << exposed(v_L) << exposed(v_R) << std::endl;
    #endif

    auto [P_L, EP_L] = find_path(s, v_L, nullptr);
    auto [P_R, EP_R] = find_path(t, v_R, nullptr);
    std::reverse(P_L.begin(), P_L.end());
    std::reverse(EP_L.begin(), EP_L.end());
    std::vector<NetworKit::node> P;
    std::vector<NetworKit::node> EP;
    for (auto v : P_L) P.push_back(v);
    for (auto e : EP_L) EP.push_back(e);
    EP.push_back(e);
    for (auto v : P_R) P.push_back(v);
    for (auto e : EP_R) EP.push_back(e);

    #if DEBUG_LOGGING
    std::cerr << "AUGMENTING ON PATH " << v_L << " ~~ " << s << " - " << t << " ~~ " << v_R << " :" << std::endl;
    print_path(P, EP);
    #endif

    for (int i = 0; i < P.size(); i += 2) {
        auto u = P[i];
        auto v = P[i + 1];
        V[u].match = v;
        V[v].match = u;
        V[u].match_edge = V[v].match_edge = EP[i];
        E[EP[i]].matched = true;
        if (i + 2 < P.size()) E[EP[i+1]].matched = false;

        #if DEBUG_LOGGING
        std::cerr << "MATCH " << u << " " << v << " (" << E[EP[i]].u << " " << E[EP[i]].v << ")" << std::endl;
        if (i + 2 < P.size()) std::cerr << "UNMATCH " << " (" << E[EP[i+1]].u << " " << E[EP[i+1]].v << ")" << std::endl;
        #endif
    }
    erase(P);
}

void MicaliVaziraniMatching::left_dfs(
    NetworKit::node s, NetworKit::node& v_L, NetworKit::node& v_R, 
    NetworKit::node& dcv, NetworKit::node& barrier) {
    
    for (auto [u, e] : V[v_L].predecessors) {
        if (E[e].used || V[u].erased) continue;
        E[e].used = true;

        #if DEBUG_LOGGING
        std::cerr << "USE EDGE (" << E[e].u << ", " << E[e].v << ")\n";
        #endif

        if (V[u].bloom != nullptr) {
            u = base_star(V[u].bloom);
        }
        if (V[u].mark == VertexData::Mark::unmarked) {
            V[u].mark = VertexData::Mark::left;
            V[u].parent = v_L;
            V[u].parent_edge = e;
            v_L = u;
            bloom_nodes.push_back(v_L);
            return;
        }
    }

    if (v_L == s) {
        bloom_found = true;
    } else {
        v_L = V[v_L].parent;
    }
}

void MicaliVaziraniMatching::right_dfs(
    NetworKit::node& v_L, NetworKit::node& v_R, 
    NetworKit::node& dcv, NetworKit::node& barrier) {

    for (auto [u, e] : V[v_R].predecessors) {
        if (E[e].used || V[u].erased) continue;
        E[e].used = true;

        #if DEBUG_LOGGING
        std::cerr << "USE EDGE (" << E[e].u << ", " << E[e].v << ")\n";
        #endif

        if (V[u].bloom != nullptr) {
            u = base_star(V[u].bloom);
        }
        if (V[u].mark == VertexData::Mark::unmarked) {
            V[u].mark = VertexData::Mark::right;
            V[u].parent = v_R;
            V[u].parent_edge = e;
            v_R = u;
            bloom_nodes.push_back(v_R);
            return;
        } else if (u == v_L) {
            dcv = u;
        }
    }

    if (v_R == barrier) {
        v_R = dcv;
        barrier = dcv;
        V[v_R].mark = VertexData::Mark::right;
        v_L = V[v_L].parent;
        bloom_nodes.push_back(v_R);
    } else {
        v_R = V[v_R].parent;
    }
}

void MicaliVaziraniMatching::erase(std::vector<NetworKit::node>& Y) {
    for (int i = 0; i < Y.size(); ++ i) {
        auto y = Y[i];
        V[y].erased = true;

        #if DEBUG_LOGGING
        std::cerr << "ERASE " << y << std::endl;
        #endif

        for (auto z : V[y].successors) {
            if (V[z].erased) continue;

            V[z].count --;
            if (V[z].count == 0) {
                Y.push_back(z);
            }
        }
    }
}

std::pair<std::list<NetworKit::node>, std::list<NetworKit::edgeid>> 
MicaliVaziraniMatching::find_path(
    NetworKit::node high, NetworKit::node low, Bloom* B) {

    #if DEBUG_LOGGING
    std::cerr << "FIND PATH " << high << " " << low << " in " << B << std::endl;
    #endif 

    if (high == low) {
        return {{ high }, {}};
    }

    auto v = high;
    auto u = high;
    NetworKit::edgeid ue;
    while (u != low) {
        #if DEBUG_LOGGING
        std::cerr << "LOOK AT " << v << std::endl;
        #endif 

        bool unvisited = false;
        for (auto [q, e] : V[v].predecessors) {
            if (E[e].visited) continue;

            #if DEBUG_LOGGING
            std::cerr << "FOUND UNVISITED PRED " << u <<  "\n";
            #endif

            if (V[q].bloom == nullptr || V[q].bloom == B) {
                E[e].visited = true;

                #if DEBUG_LOGGING
                std::cerr << "MARK (" << E[e].u << ", " << E[e].v <<  ") visited\n";
                #endif

                u = q;
                ue = e;
            } else {
                u = V[v].bloom->base;
                ue = NetworKit::none;
            }

            unvisited = true;
            break;
        }

        if (!unvisited) {
            v = V[v].parent;
        } else if (!V[u].erased && level(u) >= level(low) && 
                    (u == low || (!V[u].visited && V[u].mark == V[high].mark))) {
            
            #if DEBUG_LOGGING
            std::cerr << "PARENT(" << u << ") = " << v << std::endl;
            #endif 

            V[u].visited = true;
            V[u].parent = v;
            V[u].parent_edge = ue;
            v = u;
        }
    }

    std::list<NetworKit::node> P;
    std::list<NetworKit::edgeid> EP;
    v = low;
    do {
        #if DEBUG_LOGGING
        std::cerr << "APPEND " << v << " WITH PARENT " << V[v].parent << std::endl;
        #endif 

        P.push_back(v);
        EP.push_back(V[v].parent_edge);
        v = V[v].parent;
    } while (v != high);
    P.push_back(v);
    std::reverse(P.begin(), P.end());
    std::reverse(EP.begin(), EP.end());

    #if DEBUG_LOGGING
    std::cerr << "PREPATH:\n";
    print_path(P, EP);
    #endif 

    auto it = P.begin();
    auto eit = EP.begin();
    while (std::next(it) != P.end() && it != P.end()) {
        auto x = it;
        auto y = std::next(it);
        auto e = eit;
        if (V[*x].bloom != nullptr && V[*x].bloom != B) {
            #if DEBUG_LOGGING
            std::cerr << "OPENING " << *x << " IN " << V[*x].bloom << " NOT IN " << B << std::endl;
            #endif 

            auto [O, EO] = open(*x);

            #if DEBUG_LOGGING
            std::cerr << "OPEN RESULT:\n";
            print_path(O, EO);
            #endif

            it = std::next(y);
            eit = std::next(e);
            P.splice(x, O);

            #if DEBUG_LOGGING
            std::cerr << "AFTER SPLICE:\n";
            print_nodes(P);
            std::cerr << "LOOKING AT " << *it << " AFTER " << *std::prev(it) << std::endl;
            #endif 

            P.erase(x);
            P.erase(y);
            EP.splice(e, EO);
            EP.erase(e);

            #if DEBUG_LOGGING
            std::cerr << "PATH AFTER OPENING:\n";
            print_path(P, EP);
            std::cerr << "LOOKING AT " << *it << " AFTER " << *std::prev(it) << std::endl;
            #endif 
        } else {
            it = std::next(it);
            eit = std::next(eit);
        }
    }

    return {P, EP};
}

std::pair<std::list<NetworKit::node>, std::list<NetworKit::edgeid>>  
MicaliVaziraniMatching::open(NetworKit::node x) {
    #if DEBUG_LOGGING
    std::cerr << "OPEN " << x << std::endl;
    #endif

    auto B = V[x].bloom;
    auto b = B->base;
    if (outer(x)) {
        return find_path(x, b, B);
    } else {
        auto lp = B->left_peak;
        auto rp = B->right_peak;
        if (V[x].mark == VertexData::Mark::left) {
            auto [P, EP] = find_path(lp, x, B);
            auto [P2, EP2] = find_path(rp, b, B);
            P.splice(P.end(), P2);
            EP.push_back(B->peak_edge);
            EP.splice(EP.end(), EP2);
            return {P, EP};
        } else if (V[x].mark == VertexData::Mark::right) {
            auto [P, EP] = find_path(rp, x, B);
            auto [P2, EP2] = find_path(lp, b, B);
            P.splice(P.end(), P2);
            EP.push_back(B->peak_edge);
            EP.splice(EP.end(), EP2);
            return {P, EP};
        }
    }
}

NetworKit::node MicaliVaziraniMatching::base_star(Bloom* bloom) {
    return bloom_bases.find(bloom->id);
}

bool MicaliVaziraniMatching::exposed(NetworKit::node vertex) {
    return V[vertex].match == NetworKit::none;
}

int MicaliVaziraniMatching::level(NetworKit::node vertex) {
    return std::min(V[vertex].odd_level, V[vertex].even_level);
}

bool MicaliVaziraniMatching::outer(NetworKit::node vertex) {
    return level(vertex) % 2 == 0;
}

bool MicaliVaziraniMatching::inner(NetworKit::node vertex) {
    return level(vertex) % 2 == 1;
}

void MicaliVaziraniMatching::print_path(
    const std::list<NetworKit::node>& P, const std::list<NetworKit::edgeid>& EP) {

    auto it = P.begin();
    auto eit = EP.begin();

    while (it != P.end()) {
        std::cerr << *it;
        if (eit != EP.end()) {
            if (*eit == NetworKit::none)
                std::cerr << " ~";
            else 
                std::cerr << " (" << E[*eit].u << ", " << E[*eit].v << ")";
        }
        std::cerr << " ";

        it ++;
        if (eit != EP.end()) eit ++;
    }

    std::cerr << std::endl;
}

void MicaliVaziraniMatching::print_path(
    const std::vector<NetworKit::node>& P, const std::vector<NetworKit::edgeid>& EP) {
        
    auto it = P.begin();
    auto eit = EP.begin();

    while (it != P.end()) {
        std::cerr << *it;
        if (eit != EP.end()) {
            if (*eit == NetworKit::none)
                std::cerr << " ~";
            else 
                std::cerr << " (" << E[*eit].u << ", " << E[*eit].v << ")";
        }
        std::cerr << " ";

        it ++;
        if (eit != EP.end()) eit ++;
    }

    std::cerr << std::endl;
}

void MicaliVaziraniMatching::print_state() {
    graph.forNodes([this] (NetworKit::node vertex) {
        std::cerr << "VERTEX " << vertex << ": " << std::endl;
        std::cerr << "   MATCH       " << V[vertex].match << std::endl;
        std::cerr << "   MATCH EDGE  " << V[vertex].match_edge << std::endl;
        std::cerr << "   EVEN LEVEL  " << V[vertex].even_level << std::endl;
        std::cerr << "   ODD LEVEL   " << V[vertex].odd_level << std::endl;
        std::cerr << "   PARENT      " << V[vertex].parent << std::endl;
        std::cerr << "   PARENT EDGE " << V[vertex].parent_edge << std::endl;
        std::cerr << "   PRED        "; print_nodes(V[vertex].predecessors);
        std::cerr << "   SUCC        "; print_nodes(V[vertex].successors);
        std::cerr << "   BLOOM       " << V[vertex].bloom << std::endl;
    });
}

} /* namespace Koala */
