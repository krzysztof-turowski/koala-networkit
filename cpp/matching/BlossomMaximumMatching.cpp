#include <matching/MaximumMatching.hpp>

namespace Koala {

BlossomMaximumMatching::BlossomMaximumMatching(
    NetworKit::Graph &graph, bool perfect, InitializationStrategy initialization):
        MaximumWeightMatching(graph, perfect),
        finished(false),
        graph_edges(graph.upperEdgeIdBound()),
        node_iter(graph.upperNodeIdBound()),
        is_in_matching(graph.upperEdgeIdBound(), false),
        edges_in_matching(0),
        matched_vertex(graph.upperNodeIdBound(), NetworKit::none),
        matched_edge(graph.upperNodeIdBound(), NetworKit::none),
        trivial_blossom(graph.upperNodeIdBound(), nullptr),
        y(graph.upperNodeIdBound()) {
    graph.forEdges([this] (
            NetworKit::node u, NetworKit::node v,
            NetworKit::edgeweight weight, NetworKit::edgeid id) {
        graph_edges[id] = {u, v, 2 * static_cast<MaximumWeightMatching::weight>(weight)};
    });

    max_weight = std::numeric_limits<weight>::min();
    for (auto [u, v, w] : graph_edges)
        max_weight = std::max(w, max_weight);

    graph.forNodes([this] (NetworKit::node vertex) {
        Blossom* blossom = Blossom::trivial(vertex);
        trivial_blossom[vertex] = blossom;
        node_iter[vertex] = blossom->nodes.begin();
        add_blossom(blossom);
    });

    initialize(initialization);
}

BlossomMaximumMatching::~BlossomMaximumMatching() {
    for (auto b : blossoms) {
        b->delete_all_children();
        delete b;
    }
}

void BlossomMaximumMatching::initialize(
    BlossomMaximumMatching::InitializationStrategy initialization) {
    switch (initialization) {
    case empty:
        graph.forNodes([this] (NetworKit::node v) {
            y[v] = max_weight / 2;
        });

        return;
    case greedy:
        graph.forNodes([this] (NetworKit::node v) {
            y[v] = 0;
        });

        graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
            auto [_, __, w] = graph_edges[e];
            std::get<2>(graph_edges[e]) = 2 * w;
            y[u] = std::max(y[u], w);
            y[v] = std::max(y[v], w);
        });

        graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeid e) {
            auto [_, __, w] = graph_edges[e];

            if (y[u] + y[v] > w ||
                matched_vertex[u] != NetworKit::none || matched_vertex[v] != NetworKit::none)
                return;

            swap_edge_in_matching(e);
        });

        return;
    }
}

void BlossomMaximumMatching::run() {
    while (!finished) {
        run_stage();
    }

    // Copy root blossom pointers to separate container
    std::vector<Blossom*> final_blossoms(blossoms.begin(), blossoms.end());
    for (auto blossom : final_blossoms) {
        expand_final_blossom(blossom);
    }

    graph.forNodes([this] (NetworKit::node u) {
        matching[u] = matched_vertex[u];
    });

    hasRun = true;
}

void BlossomMaximumMatching::run_stage() {
    initialize_stage();

    // Perform searches and adjust dual variables until
    // an augmenting path has been found or we reached the solution
    while (!finished && run_substage()) {
        adjust_dual_variables();
    }

    finish_stage();

    if (perfect && is_matching_perfect()) finished = true;
}

bool BlossomMaximumMatching::run_substage() {
    initialize_substage();

    while (has_useful_edges()) {
        auto edge = get_useful_edge();

        if (consider_edge(edge))
            return false;
    }

    return true;
}

void BlossomMaximumMatching::expand_final_blossom(Blossom* blossom) {
    auto subblossoms = blossom->subblossoms;
    expand_even_blossom(blossom);
    for (auto [b, e] : subblossoms) {
        expand_final_blossom(b);
    }
}

bool BlossomMaximumMatching::consider_edge(Edge edge) {
    auto [u, v, id] = edge;
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);

    if (u_blossom == v_blossom || is_in_matching[id]) return false;

    if (u_blossom->label != even) {
        std::swap(u, v);
        std::swap(u_blossom, v_blossom);
    }

    if (v_blossom->label == free) {
        // Useful edge to a free blossom
        // Label it and the blossom matched to it's base

        auto b = v_blossom->base;
        auto c = matched_vertex[b];

        auto c_blossom = get_blossom(c);

        v_blossom->label = odd;
        v_blossom->backtrack_edge = { u, v, id };

        c_blossom->label = even;
        c_blossom->backtrack_edge = { b, c, matched_edge[b] };

        handle_grow(v_blossom, c_blossom);
    } else if (v_blossom->label == even) {
        // Useful edge between two even vertices
        // Either an augmenting path has been found, in which case we finish the stage
        // or a new blossom is created
        if (backtrack(u_blossom, v_blossom, {u, v, id})) return true;
    }

    return false;
}

bool BlossomMaximumMatching::backtrack(Blossom* u, Blossom* v, Edge edge) {
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
    augment_path(u, v, edge, u_path, v_path);

    return true;
}

bool BlossomMaximumMatching::backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path) {
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

void BlossomMaximumMatching::add_blossom(Blossom* b) {
    blossoms.push_back(b);
    b->list_it = std::prev(blossoms.end());
}

void BlossomMaximumMatching::remove_blossom(Blossom* b) {
    blossoms.erase(b->list_it);
}

void BlossomMaximumMatching::cut_path_at(
        std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut2) {
    int index = -1;
    for (unsigned i = 0; i < path.size(); ++i) {
        if (path[i].blossom == cut || path[i].blossom == cut2) {
            index = i;
            break;
        }
    }
    if (index != -1) path.resize(index + 1);
}

void BlossomMaximumMatching::create_new_blossom(
        Blossom* u, Blossom* v, Edge edge,
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    cut_path_at(u_path, v_path.size() > 0 ? v_path.back().blossom : v, v);
    cut_path_at(v_path, u_path.size() > 0 ? u_path.back().blossom : u, u);
    if (v_path.size() > 0 && v_path.back().blossom == u) u_path.clear();
    if (u_path.size() > 0 && u_path.back().blossom == v) v_path.clear();

    std::list<std::pair<Blossom*, Edge>> subblossoms;

    for (int i = v_path.size() - 1; i > 0; --i) {
        subblossoms.emplace_back(v_path[i-1].blossom, v_path[i].edge);
    }

    if (v_path.size() > 0)
        subblossoms.emplace_back(v, v_path[0].edge);
    subblossoms.emplace_back(u, reverse(edge));

    for (int i = 0; i < u_path.size(); ++i) {
        subblossoms.emplace_back(u_path[i].blossom, reverse(u_path[i].edge));
    }

    Blossom* new_blossom = Blossom::nontrivial(subblossoms);

    for (auto [b, edge] : subblossoms) {
        remove_blossom(b);
        b->parent = new_blossom;
    }
    add_blossom(new_blossom);
    handle_new_blossom(new_blossom);

    // Only update the blossom list after the handle_new_blossom call so that
    // lists for individual subblosssoms are available there
    std::list<NetworKit::node> nodes;
    for (auto [b, edge] : subblossoms) {
        nodes.splice(nodes.end(), std::move(b->nodes));
    }
    new_blossom->nodes = std::move(nodes);
}

void BlossomMaximumMatching::swap_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];

    if (is_in_matching[edge]) {
        is_in_matching[edge] = false;
        edges_in_matching--;
        if (matched_vertex[u] == v) matched_vertex[u] = matched_edge[u] = NetworKit::none;
        if (matched_vertex[v] == u) matched_vertex[v] = matched_edge[v] = NetworKit::none;
    } else {
        is_in_matching[edge] = true;
        edges_in_matching++;
        matched_vertex[u] = v;
        matched_vertex[v] = u;
        matched_edge[u] = matched_edge[v] = edge;
    }
}

void BlossomMaximumMatching::check_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];

    if (is_in_matching[edge]) {
        matched_vertex[u] = v;
        matched_vertex[v] = u;
        matched_edge[u] = matched_edge[v] = edge;
    }
}

std::list<BlossomMaximumMatching::Blossom*> BlossomMaximumMatching::blossoms_containing(
        NetworKit::node u, Blossom* until) {
    std::list<Blossom*> res;
    Blossom* iter = trivial_blossom[u];
    while (iter != until) {
        res.push_front(iter);
        iter = iter->parent;
    }
    return res;
}

void BlossomMaximumMatching::augment_path(
        Blossom* u, Blossom* v, Edge edge,
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    for (auto [b, e] : u_path) swap_edge_in_matching(e.id);
    for (auto [b, e] : v_path) swap_edge_in_matching(e.id);
    swap_edge_in_matching(edge.id);

    if (u_path.size() > 0) {
        Blossom* x = u_path[u_path.size() - 1].blossom;
        swap_edges_on_even_path(x, x->base, u_path[u_path.size() - 1].edge.u);

        for (int i = 0; i < u_path.size() - 1; ++i) {
            swap_edges_on_even_path(u_path[i].blossom, u_path[i].edge.u, u_path[i+1].edge.v);
        }
    }

    swap_edges_on_even_path(u, edge.u, u->base);
    swap_edges_on_even_path(v, edge.v, v->base);

    if (v_path.size() > 0) {
        for (int i = 0; i < v_path.size() - 1; ++i) {
            swap_edges_on_even_path(v_path[i].blossom, v_path[i].edge.u, v_path[i+1].edge.v);
        }

        Blossom* y = v_path[v_path.size() - 1].blossom;
        swap_edges_on_even_path(y, y->base, v_path[v_path.size() - 1].edge.u);
    }
}

void BlossomMaximumMatching::swap_edges_on_even_path(
        Blossom* blossom, NetworKit::node u, NetworKit::node v) {
    if (blossom->is_trivial() || u == v) return;
    auto out_vertex = u == blossom->base ? v : u;
    blossom->base = out_vertex;
    blossom->base_blossoms = blossoms_containing(out_vertex, blossom);
}

void BlossomMaximumMatching::swap_edges_on_even_path(
        Blossom* blossom, NetworKit::node out_vertex, std::list<Blossom*>&& base_blossoms) {
    if (blossom->is_trivial()) return;
    blossom->base = out_vertex;
    blossom->base_blossoms = std::move(base_blossoms);
}

void BlossomMaximumMatching::lazy_augment_path_in_blossom(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    if (blossom->base == blossom->initial_base) {
        for (auto [b, e] : blossom->subblossoms) check_edge_in_matching(e.id);
        return;
    }

    auto out_vertex = blossom->base;
    Blossom* out_blossom = blossom->base_blossoms.front();
    blossom->base_blossoms.pop_front();
    auto [pathA, pathB] = split_subblossoms(blossom->subblossoms, out_blossom);
    Blossom* base_blossom = std::get<Blossom*>(blossom->subblossoms.back());

    if (pathA.size() % 2 == 0) {
        for (auto [b, e] : pathA) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathB) check_edge_in_matching(e.id);

        if (out_blossom != base_blossom)
            swap_edges_on_even_path(base_blossom, blossom->initial_base, (pathA.begin()->second).u);

        for (auto iter = pathA.begin(); std::next(iter) != pathA.end(); iter++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(blossom->base_blossoms));
    } else {
        for (auto [b, e] : pathB) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathA) check_edge_in_matching(e.id);

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(blossom->base_blossoms));

        for (auto iter = pathB.begin(); std::next(iter) != pathB.end(); iter++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        if (out_blossom != base_blossom)
            swap_edges_on_even_path(base_blossom, blossom->initial_base, pathB.back().second.v);
    }

    handle_subblossom_shift(blossom, out_blossom);

    blossom->nodes.splice(blossom->nodes.end(), blossom->nodes,
        blossom->nodes.begin(), std::next(node_iter[pathA.back().first->last_node]));

    pathB.splice(pathB.end(), pathA);
    blossom->subblossoms = std::move(pathB);
}

std::pair<BlossomMaximumMatching::Blossom::SubblossomList,
          BlossomMaximumMatching::Blossom::SubblossomList>
BlossomMaximumMatching::split_subblossoms(Blossom::SubblossomList subblossoms, Blossom* blossom) {
    auto iter = subblossoms.begin();
    while (iter->first != blossom) iter++;
    return { { subblossoms.begin(), std::next(iter) }, { std::next(iter), subblossoms.end() } };
}

bool BlossomMaximumMatching::is_exposed(BlossomMaximumMatching::Blossom* b) {
    return matched_vertex[b->base] == NetworKit::none;
}

void BlossomMaximumMatching::adjust_dual_variables() {
    auto delta1 = perfect ? infinite_weight : calc_delta1();  // min u_i : i - even vertex
    auto delta2 = calc_delta2();  // min pi_ij : i - even vertex, j - free vertex
    auto delta3 = calc_delta3();  // min pi_ij / 2 : i,j - even vertices in different blossoms
    auto delta4 = calc_delta4();  // min z_k : B_k - odd blossom

    auto delta = std::min(delta1, std::min(delta2, std::min(delta3, delta4)));

    adjust_by_delta(delta);

    if (!perfect && delta == delta1) {
        // Optimal solution has been reached
        finished = true;
    } else if (delta == delta3) {
        find_delta3_useful_edges();
    } else if (delta == delta2) {
        find_delta2_useful_edges();
    } else if (delta == delta4) {
        // Blossoms where minimum was achieved need to be expanded
        std::vector<Blossom*> to_expand = get_odd_blossoms_to_expand();
        for (auto b : to_expand) expand_odd_blossom(b);
    }
}

void BlossomMaximumMatching::expand_odd_blossom(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    lazy_augment_path_in_blossom(blossom);

    auto node = blossom->backtrack_edge.v;
    Blossom* inner_blossom = blossoms_containing(node, blossom).front();
    Blossom* base_blossom = blossom->subblossoms.back().first;
    auto [pathA, pathB] = split_subblossoms(blossom->subblossoms, inner_blossom);

    for (auto [b, e] : blossom->subblossoms) b->label = free;

    inner_blossom->label = odd;
    inner_blossom->backtrack_edge = blossom->backtrack_edge;

    if (pathB.size() % 2 == 0) {
        int parity = 0;
        for (auto iter = pathB.begin(); iter != pathB.end(); iter++, parity++) {
            iter->first->label = parity % 2 == 0 ? even : odd;
            iter->first->backtrack_edge = iter->second;
        }
    } else {
        int parity = 0;
        Blossom* prev = base_blossom;
        for (auto iter = pathA.begin(); iter != pathA.end(); iter++, parity++) {
            prev->label = parity % 2 == 0 ? odd : even;
            prev->backtrack_edge = reverse(iter->second);
            prev = iter->first;
        }
    }

    for (auto [b, e] : blossom->subblossoms) {
        b->parent = nullptr;
        b->nodes.splice(b->nodes.end(), blossom->nodes,
            blossom->nodes.begin(), std::next(node_iter[b->last_node]));
        add_blossom(b);
    }
    remove_blossom(blossom);

    handle_odd_blossom_expansion(blossom);
    delete blossom;
}

void BlossomMaximumMatching::expand_even_blossom(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    lazy_augment_path_in_blossom(blossom);

    for (auto [b, e] : blossom->subblossoms) {
        b->parent = nullptr;
        b->nodes.splice(b->nodes.end(), blossom->nodes,
            blossom->nodes.begin(), std::next(node_iter[b->last_node]));
        add_blossom(b);
    }
    handle_even_blossom_expansion(blossom);
    remove_blossom(blossom);
    delete blossom;
}

bool BlossomMaximumMatching::is_matching_perfect() {
    return 2 * edges_in_matching == graph.upperNodeIdBound();
}

void BlossomMaximumMatching::Blossom::delete_all_children() {
    for (auto [b, e] : subblossoms) {
        b->delete_all_children();
        delete b;
    }
}

BlossomMaximumMatching::Blossom::~Blossom() {
    if (data != nullptr) delete data;
}

BlossomMaximumMatching::Blossom*
BlossomMaximumMatching::Blossom::trivial(NetworKit::node vertex) {
    return new Blossom {
        nullptr, vertex, vertex, vertex, {}, {}, { vertex },
        free, no_edge, false,
        0, {}, nullptr
    };
}

BlossomMaximumMatching::Blossom*
BlossomMaximumMatching::Blossom::nontrivial(const SubblossomList& subblossoms) {
    Blossom* base = subblossoms.back().first;
    return new Blossom {
        nullptr, base->base, base->base, base->last_node, {}, subblossoms, {},
        even, base->backtrack_edge, false, 0, {}, nullptr
    };
}

bool BlossomMaximumMatching::Blossom::is_trivial() {
    return subblossoms.size() == 0;
}

void BlossomMaximumMatching::Blossom::for_nodes(
        const std::function<void(NetworKit::node)>& handle) {
    if (is_trivial()) {
        handle(base);
    } else {
        for (auto [b, edge] : subblossoms)
            b->for_nodes(handle);
    }
}

bool BlossomMaximumMatching::Blossom::contains(NetworKit::node v) {
    if (is_trivial()) return base == v;

    for (auto [b, e] : subblossoms) {
        if (b->contains(v)) return true;
    }
    return false;
}

BlossomMaximumMatching::Edge BlossomMaximumMatching::reverse(const Edge& info) {
    auto [u, v, id] = info;
    return { v, u, id };
}

} /* namespace Koala */
