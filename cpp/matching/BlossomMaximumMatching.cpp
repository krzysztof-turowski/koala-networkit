#include <matching/MaximumMatching.hpp>

namespace Koala {

BlossomMaximumMatching::BlossomMaximumMatching(NetworKit::Graph &graph):
        MaximumWeightMatching(graph),
        finished(false),
        graph_edges(graph.upperEdgeIdBound(), {NetworKit::none, NetworKit::none, 0}),
        node_iter(graph.upperNodeIdBound()),
        is_in_matching(graph.upperEdgeIdBound(), false),
        matched_vertex(graph.upperNodeIdBound(), NetworKit::none),
        matched_edge(graph.upperNodeIdBound(), NetworKit::none),
        trivial_blossom(graph.upperNodeIdBound(), nullptr) {
    
    graph.forEdges([this] (
            NetworKit::node u, NetworKit::node v,
            NetworKit::edgeweight weight, NetworKit::edgeid id) {
        graph_edges[id] = {u, v, 2 * static_cast<MaximumWeightMatching::weight>(weight)};
    });

    max_weight = std::numeric_limits<weight>::min();
    for (auto [u, v, w] : graph_edges) 
        max_weight = std::max(w, max_weight);

    graph.forNodes([this] (NetworKit::node vertex) {
        Blossom* blossom = new Blossom {
            nullptr, vertex, vertex, vertex, {}, {}, { vertex },
            free, no_edge, false,
            0, {}, nullptr
        };
        trivial_blossom[vertex] = blossom;
        node_iter[vertex] = blossom->nodes.begin();
        add_blossom(blossom);
    });
}

BlossomMaximumMatching::~BlossomMaximumMatching() {
    for (auto b : blossoms) {
        b->delete_all_children();
        delete b;
    }
}

void BlossomMaximumMatching::run() {
    while (!finished) {
        run_stage();
    }

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
    #if DEBUG_LOGGING
    std::cerr << "\n========================== START  STAGE ==========================\n";
    #endif

    initialize_stage();

    #if DEBUG_LOGGING
    check_consistency();
    #endif

    // Perform searches and adjust dual variables until 
    // an augmenting path has been found or we reached the solution
    while (!finished && run_substage()) {
        adjust_dual_variables();

        #if DEBUG_LOGGING
        check_consistency();
        #endif
    }

    #if DEBUG_LOGGING
    check_consistency();
    #endif

    finish_stage();
}

bool BlossomMaximumMatching::run_substage() {
    #if DEBUG_LOGGING
    std::cerr << "\n========================== START SUBSTAGE ==========================\n";
    #endif

    initialize_substage();

    while (has_useful_edges()) {
        auto [u, v, id] = get_useful_edge();

        if (consider_edge({u, v, id})) return false;

        #if DEBUG_LOGGING
        check_consistency();
        #endif
    }

    return true;
}

void BlossomMaximumMatching::expand_final_blossom(Blossom* blossom) {
    #if DEBUG_LOGGING
    std::cerr << "Expand final blossom "; blossom->nodes_print(); std::cerr << std::endl;
    #endif

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

    #if DEBUG_LOGGING
    std::cerr << "> Consider edge (" << u << ", " << v << ")\n";
    #endif

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

        #if DEBUG_LOGGING
        std::cerr << "STEP GROW\n";
        v_blossom->nodes_print();
        std::cerr << ", "; c_blossom->nodes_print();
        std::cerr << std::endl;
        #endif

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
    for (unsigned i = 0; i < path.size(); ++ i) {
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

    #if DEBUG_LOGGING
    std::cerr << "STEP BLOSSOM" << std::endl;
    print_backtrack(u, v, edge, u_path, v_path);
    #endif

    Blossom* base = v_path.size() > 0 ? v_path.back().blossom : v;

    std::list<std::pair<Blossom*, Edge>> subblossoms;
    
    for (int i = v_path.size() - 1; i > 0; -- i) {
        subblossoms.emplace_back(v_path[i-1].blossom, v_path[i].edge);
    }

    if (v_path.size() > 0) 
        subblossoms.emplace_back(v, v_path[0].edge);
    subblossoms.emplace_back(u, reverse(edge));

    for (int i = 0; i < u_path.size(); ++ i) {
        subblossoms.emplace_back(u_path[i].blossom, reverse(u_path[i].edge));
    }

    Blossom* new_blossom = new Blossom {
        nullptr, base->base, base->base, base->last_node, {}, subblossoms, {},
        even, base->backtrack_edge, false, 0, {}, nullptr
    };

    for (auto [b, edge] : subblossoms) { 
        remove_blossom(b);
        b->parent = new_blossom;
    }
    add_blossom(new_blossom);
    handle_new_blossom(new_blossom);

    std::list<NetworKit::node> nodes;
    for (auto [b, edge] : subblossoms) { 
        nodes.splice(nodes.end(), std::move(b->nodes));
    }
    new_blossom->nodes = std::move(nodes);

    #if DEBUG_LOGGING
    std::cerr << "Created new blossom: " << std::endl;
    new_blossom->short_print(); std::cerr << std::endl;
    #endif
}

void BlossomMaximumMatching::swap_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    
    if (is_in_matching[edge]) {
        #if DEBUG_LOGGING
        std::cerr << "Remove edge (" << u << ", " << v << ")" << std::endl;
        #endif

        is_in_matching[edge] = false;
        if (matched_vertex[u] == v) matched_vertex[u] = matched_edge[u] = NetworKit::none;
        if (matched_vertex[v] == u) matched_vertex[v] = matched_edge[v] = NetworKit::none;
    } else {
        #if DEBUG_LOGGING
        std::cerr << "Add edge (" << u << ", " << v << ")" << std::endl;
        #endif

        is_in_matching[edge] = true;
        matched_vertex[u] = v;
        matched_vertex[v] = u;
        matched_edge[u] = matched_edge[v] = edge;
    }
}

void BlossomMaximumMatching::check_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    
    if (is_in_matching[edge]) {
        #if DEBUG_LOGGING
        std::cerr << "Ensure edge (" << u << ", " << v << ")" << std::endl;
        #endif

        matched_vertex[u] = v;
        matched_vertex[v] = u;
        matched_edge[u] = matched_edge[v] = edge;
    }
}

std::list<BlossomMaximumMatching::Blossom*> BlossomMaximumMatching::blossoms_containing(
        NetworKit::node u, Blossom* until) {
    #if DEBUG_LOGGING
    std::cerr << "find blossoms on path to " << u << " until blossom ";
    until->nodes_print(); std::cerr << std::endl;
    #endif

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
    #if DEBUG_LOGGING
    std::cerr << "STEP AUGMENT" << std::endl;
    print_backtrack(u, v, edge, u_path, v_path);
    #endif
    
    for (auto [b, e] : u_path) swap_edge_in_matching(e.id);
    for (auto [b, e] : v_path) swap_edge_in_matching(e.id);
    swap_edge_in_matching(edge.id);

    if (u_path.size() > 0) {
        Blossom* x = u_path[u_path.size() - 1].blossom;
        swap_edges_on_even_path(x, x->base, u_path[u_path.size() - 1].edge.u);

        for (int i = 0; i < u_path.size() - 1; ++ i) {
            swap_edges_on_even_path(u_path[i].blossom, u_path[i].edge.u, u_path[i+1].edge.v);
        }
    }

    swap_edges_on_even_path(u, edge.u, u->base);
    swap_edges_on_even_path(v, edge.v, v->base);

    if (v_path.size() > 0) {
        for (int i = 0; i < v_path.size() - 1; ++ i) {
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

    #if DEBUG_LOGGING
    std::cerr << "Schedule swap of edges on the even path " << blossom->initial_base << " ~~> " << out_vertex
              << " in blossom "; blossom->short_print(); std::cerr << std::endl;
    #endif
}

void BlossomMaximumMatching::swap_edges_on_even_path(
        Blossom* blossom, NetworKit::node out_vertex, std::list<Blossom*>&& base_blossoms) {
    
    if (blossom->is_trivial()) return;
    blossom->base = out_vertex;
    blossom->base_blossoms = std::move(base_blossoms);
   
    #if DEBUG_LOGGING
    std::cerr << "Schedule deeper swap of edges on the even path " << blossom->base << " ~~> " 
              << out_vertex << " in blossom "; blossom->nodes_print(); std::cerr << std::endl;
    #endif
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

    #if DEBUG_LOGGING
    std::cerr << "Lazy augment " << out_vertex << " inside "; out_blossom->nodes_print(); 
    std::cerr << "\nin blossom "; blossom->nodes_print(); std::cerr<<std::endl;
    #endif

    if (pathA.size() % 2 == 0) {
        for (auto [b, e] : pathA) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathB) check_edge_in_matching(e.id);

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, blossom->initial_base, (pathA.begin()->second).u);

        for (auto iter = pathA.begin(); std::next(iter) != pathA.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(blossom->base_blossoms));
    } else {
        for (auto [b, e] : pathB) swap_edge_in_matching(e.id);
        for (auto [b, e] : pathA) check_edge_in_matching(e.id);

        swap_edges_on_even_path(out_blossom, out_vertex, std::move(blossom->base_blossoms));

        for (auto iter = pathB.begin(); std::next(iter) != pathB.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, blossom->initial_base, pathB.back().second.v);
    }

    handle_subblossom_shift(blossom, out_blossom);

    blossom->nodes.splice(blossom->nodes.end(), blossom->nodes, 
        blossom->nodes.begin(), std::next(node_iter[pathA.back().first->last_node]));
    
    pathB.splice(pathB.end(), pathA);
    blossom->subblossoms = std::move(pathB); // pathB + pathA
}

std::pair<std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::Edge>>, 
          std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::Edge>>>
BlossomMaximumMatching::split_subblossoms(
        std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::Edge>> subblossoms, 
        Blossom* blossom) {
    
    auto iter = subblossoms.begin();
    while (iter->first != blossom) iter ++;
    return { { subblossoms.begin(), std::next(iter) }, { std::next(iter), subblossoms.end() } };
}

bool BlossomMaximumMatching::is_exposed(BlossomMaximumMatching::Blossom* b) {
    return matched_vertex[b->base] == NetworKit::none;
}

void BlossomMaximumMatching::adjust_dual_variables() {
    #if DEBUG_LOGGING
    std::cerr << "STEP DUAL ADJUSTMENT\n";
    check_consistency();
    for (auto b : blossoms) b->check_nodes_list();
    #endif

    auto delta1 = calc_delta1(); // min u_i : i - even vertex
    auto delta2 = calc_delta2(); // min pi_ij : i - even vertex, j - free vertex
    auto delta3 = calc_delta3(); // min pi_ij / 2 : i,j - even vertices in different blossoms
    auto delta4 = calc_delta4(); // min z_k : B_k - odd blossom 

    auto delta = std::min(delta1, std::min(delta2, std::min(delta3, delta4)));

    #if DEBUG_LOGGING
    std::cerr << "-------- CALCULATING BEST DELTA --------\n";
    std::cerr << "delta_1 = " << delta1 << std::endl;
    std::cerr << "delta_2 = " << delta2 << std::endl;
    std::cerr << "delta_3 = " << delta3 << std::endl;
    std::cerr << "delta_4 = " << delta4 << std::endl;
    std::cerr << "Adjusting dual variables by " << delta << std::endl;
    #endif

    adjust_by_delta(delta);

    if (delta == delta1) {
        // Optimal solution has been reached
        finished = true;
    } 
    else if (delta == delta3) find_delta3_useful_edges();
    else if (delta == delta2) find_delta2_useful_edges();
    else if (delta == delta4) {
        // Blossoms where minimum was achieved need to be expanded
        std::vector<Blossom*> to_expand = get_odd_blossoms_to_expand();
        for (auto b : to_expand) expand_odd_blossom(b);
    }
}

void BlossomMaximumMatching::expand_odd_blossom(Blossom* blossom) {
    if (blossom->is_trivial()) return;

    #if DEBUG_LOGGING
    std::cerr << "Expanding odd blossom "; blossom->nodes_print(); std::cerr << std::endl;
    #endif

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
    #if DEBUG_LOGGING
    std::cerr << "Expanding even blossom "; blossom->nodes_print(); std::cerr << std::endl;
    #endif

    lazy_augment_path_in_blossom(blossom);

    if (blossom->is_trivial()) return;
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

void BlossomMaximumMatching::print_backtrack(
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

void BlossomMaximumMatching::Blossom::delete_all_children() {
    for (auto [b, e] : subblossoms) {
        b->delete_all_children();
        delete b;
    }
}

BlossomMaximumMatching::Blossom::~Blossom() {
    if (data != nullptr) delete data;
}

bool BlossomMaximumMatching::Blossom::is_trivial() {
    return subblossoms.size() == 0;
}

void BlossomMaximumMatching::Blossom::for_nodes(const std::function<void(NetworKit::node)>& handle) {
    if (is_trivial()) {
        handle(base);
    } else {
        for (auto [b, edge] : subblossoms)
            b->for_nodes(handle);
    }
}

void BlossomMaximumMatching::Blossom::check_consistency() {
    for (auto [b, e] : subblossoms) {
        auto [u, v, id] = e;
        if (!b->contains(v)) {
            std::cerr << "Subblossom edge in wrong direction "; 
            short_print(); 
            std::cerr << std::endl;
            exit(1);
        }
        b->check_consistency();
    }
}

void BlossomMaximumMatching::Blossom::check_nodes_list() {
    std::list<NetworKit::node> deep_list;
    for_nodes([&deep_list] (NetworKit::node v) {
        deep_list.push_back(v);
    });
    
    if (deep_list.size() != nodes.size()) {
        std::cerr << "Blossom list with wrong length for ";
        nodes_print();
        std::cerr << std::endl;
        for (auto n : nodes) std::cerr << n << " ";
        std::cerr << std::endl;
        exit(1);
    }

    auto it = deep_list.begin();
    for (auto n : nodes) {
        if (n != *it) {
            std::cerr << "Blossom list wrong for ";
            nodes_print();
            std::cerr << std::endl;
            for (auto n : nodes) std::cerr << n << " ";
            std::cerr << std::endl;
            exit(1);
        }
        it ++;
    }
}

bool BlossomMaximumMatching::Blossom::contains(NetworKit::node v) {
    if (is_trivial()) return base == v;

    for (auto [b, e] : subblossoms) {
        if (b->contains(v)) return true;
    }
    return false;
}

void BlossomMaximumMatching::Blossom::print(int depth) {
    std::string padding;
    for (int i = 0; i < depth; ++ i) padding += "    ";
    std::cerr << padding << "Blossom:" << std::endl;
    if (depth == 0) {
        std::cerr << padding << "    label:     ";
        std::cerr << (label == even ? "S" : label == odd ? "T" : "free") << std::endl;
    }
    std::cerr 
        << padding << "    base:       " << base << std::endl
        << padding << "    weight:     " << z << std::endl
        << padding << "    label:      " 
            << (label == even ? "even" : label == odd ? "odd" : "free") << std::endl
        << padding << "    backtrack:  " 
            << "("  << backtrack_edge.u 
            << ", " << backtrack_edge.v
            << ") id: "  << backtrack_edge.id << std::endl
        << padding << "    subblossoms:" << std::endl;
    for (auto [b, edge] : subblossoms) {
        std::cerr 
            << padding << "  >>edge:      " 
                << "("  << edge.u << ", " << edge.v << ") id: "  << edge.id << std::endl;
        b->print(depth + 1);
    }
}

void BlossomMaximumMatching::Blossom::short_print() {
    std::cerr << "{ "; 
    if (is_trivial()) {
        std::cerr << base << " ";
    } else {
        for (auto [b, e] : subblossoms) {
            std::cerr << "(";
            b->short_print();
            std::cerr << ", (" << e.u << ", " << e.v << ")) ";
        }
    }
    if (parent == nullptr) {
        std::cerr << "| " << (label == even ? "S" : label == odd ? "T" : "-") << " ";
        if (backtrack_edge.id != NetworKit::none) {
            std::cerr << "(" << backtrack_edge.u << ", " << backtrack_edge.v << ") ";
        }
    }
    std::cerr << "}";
}

void BlossomMaximumMatching::Blossom::nodes_print() {
    std::cerr << "("; 
    if (is_trivial()) {
        std::cerr << base;
    } else {
        for (auto sb : subblossoms) sb.first->nodes_print();
    }
    if (parent == nullptr) {
        std::cerr << "|" << (label == even ? "S" : label == odd ? "T" : "-") 
                  << "|" << base;
    }
    std::cerr << ")";
}

BlossomMaximumMatching::Edge BlossomMaximumMatching::reverse(const Edge& info) {
    auto [u, v, id] = info;
    return { v, u, id };
}

std::string BlossomMaximumMatching::edge_to_string(const BlossomMaximumMatching::Edge& e) { 
    return (e == no_edge) ? 
        "none" : 
        "(" + std::to_string(e.u) + ", " + std::to_string(e.v) + ")"; 
}

} /* namespace Koala */
