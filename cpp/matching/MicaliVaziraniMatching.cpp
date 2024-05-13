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
        candidates(graph.numberOfNodes()),
        bridges(2 * graph.numberOfNodes() + 1),
        bloom_bases(graph.upperNodeIdBound()) {
            graph.forNodes([this] (NetworKit::node vertex) {
                V[vertex].match = NetworKit::none;
            });
        }

MicaliVaziraniMatching::MicaliVaziraniMatching(
    NetworKit::Graph &graph, const std::vector<NetworKit::node>& initial_matching):
        MaximumCardinalityMatching(graph),
        V(graph.upperNodeIdBound()),
        E(graph.upperEdgeIdBound()),
        candidates(graph.numberOfNodes()),
        bridges(2 * graph.numberOfNodes() + 1),
        bloom_bases(graph.upperNodeIdBound()) {
            graph.forNodes([this, &initial_matching] (NetworKit::node vertex) {
                V[vertex].match = initial_matching[vertex];
            });
        }

void MicaliVaziraniMatching::reset() {
    augmentation_happened = false;
    color_counter = 1;

    graph.forNodes([this] (NetworKit::node vertex) {
        V[vertex].parent = NetworKit::none;
        V[vertex].even_level = inf_level;
        V[vertex].odd_level = inf_level;
        V[vertex].bloom = nullptr;
        V[vertex].predecessors.clear();
        V[vertex].pred_it = 0;
        V[vertex].successors.clear();
        V[vertex].children.clear();
        V[vertex].pred_count = 0;
        V[vertex].erased = false;
        V[vertex].color = no_color;
    });

    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        E[e].type = EdgeData::Type::none;
    });

    for (NetworKit::index i = 0; i < graph.numberOfNodes(); ++i) {
        candidates[i].clear();
        bridges[2 * i + 1].clear();
        bloom_bases.reset(i, i);
    }

    erase_queue.clear();
}

void MicaliVaziraniMatching::clear_blooms() {
    for (auto B : current_blooms) delete B;
    current_blooms.clear();
}

void MicaliVaziraniMatching::run() {
    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
        E[e].u = u;
        E[e].v = v;
    });

    do {
        reset();

        search();
        
        clear_blooms();
    } while (augmentation_happened);

    graph.forNodes([this] (NetworKit::node vertex) {
        matching[vertex] = V[vertex].match;
    });

    hasRun = true;
}

void MicaliVaziraniMatching::search() {
    graph.forNodes([this] (NetworKit::node vertex) {
        if (exposed(vertex)) {
            set_level(vertex, 0);
        }
    });

    iter = max_iter = 0;
    augmentation_happened = false;
    for (iter = 0; iter < candidates.size() && !augmentation_happened; iter++) {
        for (auto v : candidates[iter]) {
            graph.forEdgesOf(v, [this] (NetworKit::node v, NetworKit::node u, NetworKit::edgeid e) {
                if (V[u].erased || E[e].type != EdgeData::Type::none ||
                    (iter % 2 == 0 && u == V[v].match) || (iter % 2 == 1 && u != V[v].match))
                    return;

                if (min_level(u) >= iter + 1) {
                    E[e].type = EdgeData::Type::prop;
                    if (min_level(u) > iter + 1) {
                        set_level(u, iter + 1);
                    }
                    V[u].predecessors.push_back(v);
                    V[u].pred_count++;
                    V[v].successors.push_back(u);
                } else {
                    E[e].type = EdgeData::Type::bridge;
                    auto t = tenacity(u, v);
                    if (t < inf_level) {
                        bridges[t].push_back(e);
                    }
                }
            });
        }

        for (int i = 0; i < bridges[2 * iter + 1].size(); ++i) {
            auto e = bridges[2 * iter + 1][i];
            auto s = E[e].v;
            auto t = E[e].u;
            if (!V[s].erased && !V[t].erased)
                bloss_aug(s, t);
        }
    }
}

void MicaliVaziraniMatching::flip_edge(NetworKit::node u, NetworKit::node v) {
    if (V[u].erased || V[v].erased || V[u].match == v) return;

    V[u].match = v;
    V[v].match = u;

    V[u].erased = true;
    V[v].erased = true;
    erase_queue.push_back(u);
    erase_queue.push_back(v);
}

bool MicaliVaziraniMatching::open(NetworKit::node x, NetworKit::node x_base, NetworKit::node y) {
    if (x_base == y) {
        find_path(x, x_base);
        return true;
    }

    for (auto [pred, bloom_exit] : V[x_base].children) {
        if (bloom_exit == y || V[bloom_exit].color == V[x_base].color) {
            if (open(pred, bloom_exit, y)) {
                find_path(x, x_base);
                flip_edge(x_base, pred);
                return true;
            }
        }
    }

    return false;
}

void MicaliVaziraniMatching::find_path(NetworKit::node x, NetworKit::node y) {
    if (x == y) return;

    if (outer(x)) {
        int z = V[x].predecessors.front();

        auto w = *std::find_if(V[z].predecessors.begin(), V[z].predecessors.end(),
            [this, z] (NetworKit::node pred) {
                return base_star(pred) == base_star(z);
        });

        flip_edge(z, w);
        find_path(w, y);
    } else {
        auto [u, v, u_0, v_0] = get_bridge(x);
        auto x_base = base(x);

        open(u, u_0, x);
        flip_edge(u, v);
        open(v, v_0, x_base);
        find_path(x_base, y);
    }
}

void MicaliVaziraniMatching::bloss_aug(NetworKit::node g, NetworKit::node r) {
    if (V[g].bloom == V[r].bloom && V[g].bloom != nullptr) return;

    auto v_G = base_star(g);
    auto v_R = base_star(r);

    if (v_G == v_R) return;

    // Different colors for each bridge
    int green_color = ++color_counter;
    int red_color = ++color_counter;

    // Remember roots of dfs trees
    auto green_root = v_G;
    auto red_root = v_R;

    V[v_G].color = green_color;
    V[v_G].parent = v_G;
    V[v_R].color = red_color;
    V[v_R].parent = r;

    // Prevents green from backtracking too far
    auto barrier = v_G;

    bridge_support = { v_G, v_R };
    bloom_found = false;

    // Advance higher dfs until an augmenting path or a bloom is found
    while ((!exposed(v_G) || !exposed(v_R)) && !bloom_found) {
        if (min_level(v_R) >= min_level(v_G)) {
            red_dfs_step(v_R, red_color, v_G, r, barrier);
        } else {
            green_dfs_step(v_G, green_color, v_R, barrier);
        }
    }

    if (bloom_found) {
        // Bloom found

        Bloom *B = new Bloom { v_G, green_color, red_color, g, green_root, r, red_root };
        current_blooms.push_back(B);

        // Base not part of bloom
        V[B->base].color = no_color;

        for (auto y : bridge_support) {
            if (V[y].color == no_color) continue;

            // Assign to bloom
            V[y].bloom = B;
            bloom_bases.link(B->base, y, B->base);

            // Assign other level
            set_level(y, 2 * iter + 1 - min_level(y));

            if (inner(y)) {
                // Find possible bridges

                graph.forEdgesOf(y,
                        [this] (NetworKit::node y, NetworKit::node u, NetworKit::edgeid e) {
                    if (E[e].type != EdgeData::Type::bridge) return;

                    int t = tenacity(y, u);
                    if (t < inf_level)
                        bridges[t].push_back(e);
                });
            }
        }
        return;
    }

    // Augmenting path found
    augmentation_happened = true;

    // Find path in red tree
    open(r, red_root, v_R);
    // Flip matching on the bridge
    flip_edge(r, g);
    // Find path in green tree
    open(g, green_root, v_G);

    erase(erase_queue);
}

void MicaliVaziraniMatching::red_dfs_step(
        NetworKit::node& v_R, int red_color, NetworKit::node& v_G,
        NetworKit::node r, NetworKit::node& barrier) {
    while (V[v_R].pred_it < V[v_R].predecessors.size()) {
        auto u = V[v_R].predecessors[V[v_R].pred_it++];
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
        NetworKit::node& v_G, int green_color, NetworKit::node& v_R,
        NetworKit::node& barrier) {
    while (V[v_G].pred_it < V[v_G].predecessors.size()) {
        auto u = V[v_G].predecessors[V[v_G].pred_it++];
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

void MicaliVaziraniMatching::erase(std::vector<NetworKit::node>& Y) {
    for (size_t i = 0; i < Y.size(); ++i) {
        auto y = Y[i];

        for (auto z : V[y].successors) {
            if (V[z].erased) continue;

            V[z].pred_count--;
            if (V[z].pred_count == 0) {
                V[z].erased = true;
                Y.push_back(z);
            }
        }
    }
}

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
    return V[u].match == v ? V[u].odd_level  + V[v].odd_level  + 1
                           : V[u].even_level + V[v].even_level + 1;
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
    return V[vertex].color == B->red_color ?
        std::make_tuple(B->red_peak, B->green_peak, B->red_root, B->green_root) :
        std::make_tuple(B->green_peak, B->red_peak, B->green_root, B->red_root);
}

} /* namespace Koala */
