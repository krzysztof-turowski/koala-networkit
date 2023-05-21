#include <matching/MaximumMatching.hpp>

namespace Koala {

MaximumMatching::MaximumMatching(NetworKit::Graph &graph) : graph(graph) { }

const std::map<NetworKit::node, NetworKit::node>& MaximumMatching::getMatching() const {
    assureFinished();
    return matching;
}

BlossomMaximumMatching::BlossomMaximumMatching(NetworKit::Graph &graph):
        MaximumMatching(graph),
        is_in_matching(graph.upperEdgeIdBound(), false),
        graph_edges(graph.upperEdgeIdBound(), {NetworKit::none, NetworKit::none, 0.0}),
        matched_vertex(graph.upperNodeIdBound(), NetworKit::none),
        trivial_blossom(graph.upperNodeIdBound(), nullptr) {
    
    graph.forEdges([this] (NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight, NetworKit::edgeid id) {
        graph_edges[id] = {u, v, weight};
    });

    graph.forNodes([this] (NetworKit::node vertex) {
            Blossom* trivial_blossom = new Blossom {
            nullptr, vertex, std::list<std::pair<Blossom*, EdgeInfo>>(), 
            free, no_edge, false,
            0, nullptr
        };
        this->trivial_blossom[vertex] = trivial_blossom;
        this->blossoms.insert(trivial_blossom);
    });
}

BlossomMaximumMatching::~BlossomMaximumMatching() {
    for (auto b : blossoms) {
        b->delete_all_children();
        delete b;
    }
}

void BlossomMaximumMatching::run() {
    while (!hasRun) {
        run_stage();
    }

    graph.forNodes([this] (NetworKit::index v) {
        matching[v] = matched_vertex[v];
    });
}

void BlossomMaximumMatching::run_stage() {
    std::cerr << "\n========================== START  STAGE ==========================\n";

    initialize_stage();

    check_consistency();

    // Perform searches and adjust dual variables until 
    // an augmenting path has been found or we reached the solution
    while (!hasRun && run_substage()) {
        adjust_dual_variables();

        check_consistency();
    }

    check_consistency();

    finish_stage();
}

bool BlossomMaximumMatching::run_substage() {
    std::cerr << "\n========================== START SUBSTAGE ==========================\n";

    initialize_substage();

    while (has_useful_edges()) {
        auto [u, v, id] = get_useful_edge();

        if (consider_edge({u, v, id})) return false;
    }

    return true;
}

bool BlossomMaximumMatching::consider_edge(EdgeInfo edge) {
    auto [u, v, id] = edge;
    auto u_blossom = get_blossom(u);
    auto v_blossom = get_blossom(v);

    if (u_blossom == v_blossom || is_in_matching[id]) return false;

    if (u_blossom->label != even) {
        std::swap(u, v);
        std::swap(u_blossom, v_blossom);
    }

    std::cerr << ">>>>>>>>>>>>>> Consider edge (" << u << ", " << v << ")\n";

    if (v_blossom->label == free) {
        // Useful edge to a free blossom
        // Label it and the blossom matched to it's base

        auto b = v_blossom->base;
        auto c = matched_vertex[b];

        auto c_blossom = get_blossom(c);

        v_blossom->label = odd;
        v_blossom->backtrack_edge = { u, v, id };
        label_odd(v_blossom);

        c_blossom->label = even;
        // TODO store edgeid of matched edges to do this in O(1)
        c_blossom->backtrack_edge = { b, c, graph.edgeId(b, c) }; 
        
        label_even(c_blossom);

        // TODO move label_odd/even to one function handle_grow

    } else if (v_blossom->label == even) {
        // Useful edge between two even vertices
        // Either an augmenting path has been found, in which case we finish the stage
        // or a new blossom is created
        if (backtrack(u_blossom, v_blossom, {u, v, id})) return true;
    }

    return false;
}

bool BlossomMaximumMatching::backtrack(Blossom* u, Blossom* v, EdgeInfo edge) {
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

    // std::cerr << "backtrack " << a << " " << b << " " << iter->visited << " ";
    // iter->short_print();

    if (iter->visited) return true;

    iter->visited = true;
    return false;
}

void BlossomMaximumMatching::cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut2) {
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
        Blossom* u, Blossom* v, EdgeInfo edge, 
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {
    
    cut_path_at(u_path, v_path.size() > 0 ? v_path.rbegin()->blossom : v, v);
    cut_path_at(v_path, u_path.size() > 0 ? u_path.rbegin()->blossom : u, u);
    if (v_path.size() > 0 && v_path.rbegin()->blossom == u) u_path.clear();
    if (u_path.size() > 0 && u_path.rbegin()->blossom == v) v_path.clear();

    std::cerr << "Creating new blossom" << std::endl;
    print_backtrack(u, v, edge, u_path, v_path);

    Blossom* base = v_path.size() > 0 ? v_path.rbegin()->blossom : v;

    std::list<std::pair<Blossom*, EdgeInfo>> sub_blossoms;
    
    for (int i = v_path.size() - 1; i > 0; -- i) {
        sub_blossoms.emplace_back(v_path[i-1].blossom, v_path[i].edge);
    }

    if (v_path.size() > 0) 
        sub_blossoms.emplace_back(v, v_path[0].edge);
    sub_blossoms.emplace_back(u, reverse(edge));

    for (int i = 0; i < u_path.size(); ++ i) {
        sub_blossoms.emplace_back(u_path[i].blossom, reverse(u_path[i].edge));
    }

    Blossom* new_blossom = new Blossom {
        nullptr, base->base, sub_blossoms,
        even, base->backtrack_edge, false, 0, nullptr
    };

    for (auto [b, edge] : sub_blossoms) { 
        blossoms.erase(b);
        b->parent = new_blossom;
    }
    blossoms.insert(new_blossom);
    handle_new_blossom(new_blossom);

    std::cerr << "Created new blossom: " << std::endl;
    new_blossom->short_print(); std::cerr << std::endl;
}

void BlossomMaximumMatching::augment_path(
        Blossom* u, Blossom* v, EdgeInfo edge,
        std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path) {

    std::cerr << "Augmenting the matching on path:" << std::endl;
    print_backtrack(u, v, edge, u_path, v_path);
    
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

void BlossomMaximumMatching::swap_edge_in_matching(NetworKit::edgeid edge) {
    auto [u, v, w] = graph_edges[edge];
    
    if (is_in_matching[edge]) {
        // std::cerr << "Remove edge (" << u << ", " << v << ")" << std::endl;
        is_in_matching[edge] = false;
        if (matched_vertex[u] == v) matched_vertex[u] = NetworKit::none;
        if (matched_vertex[v] == u) matched_vertex[v] = NetworKit::none;
    } else {
        // std::cerr << "Add edge (" << u << ", " << v << ")" << std::endl;
        is_in_matching[edge] = true;
        matched_vertex[u] = v;
        matched_vertex[v] = u;
    }
}

void BlossomMaximumMatching::swap_edges_on_even_path(
        Blossom* blossom, NetworKit::node u, NetworKit::node v) {
    
    if (blossom->is_trivial() || u == v) return;
    auto out_vertex = u == blossom->base ? v : u;
    auto blossoms = blossoms_containing(out_vertex, blossom);
    auto iter = blossoms.rbegin();
    // std::cerr << "Swap edges on the even path " << blossom->base << " ~~> " << out_vertex << " in blossom "; blossom->short_print(); std::cerr << std::endl;
    augment_path_on_blossoms(blossom, iter, out_vertex);
}

void BlossomMaximumMatching::augment_path_on_blossoms(
        Blossom* blossom, std::vector<Blossom*>::reverse_iterator out_iter, NetworKit::node out_vertex) {

    if (blossom->is_trivial()) return;
    
    // std::cerr << "Augment " << out_vertex << " inside "; (*out_iter)->short_print(); std::cerr<<std::endl;
    // std::cerr << "in blossom "; blossom->short_print(); std::cerr<<std::endl;

    Blossom* out_blossom = *out_iter;
    auto [pathA, pathB] = split_subblossoms(blossom->sub_blossoms, out_blossom);
    Blossom* base_blossom = std::get<Blossom*>(*blossom->sub_blossoms.rbegin());

    if (pathA.size() % 2 == 0) {
        for (auto [b, e] : pathA) swap_edge_in_matching(e.id);

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, blossom->base, (pathA.begin()->second).u);

        for (auto iter = pathA.begin(); std::next(iter) != pathA.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        augment_path_on_blossoms(out_blossom, std::next(out_iter), out_vertex);
    } else {
        for (auto [b, e] : pathB) swap_edge_in_matching(e.id);

        augment_path_on_blossoms(out_blossom, std::next(out_iter), out_vertex);

        for (auto iter = pathB.begin(); std::next(iter) != pathB.end(); iter ++) {
            swap_edges_on_even_path(iter->first, iter->second.v, std::next(iter)->second.u);
        }

        if (out_blossom != base_blossom) 
            swap_edges_on_even_path(base_blossom, blossom->base, pathB.rbegin()->second.v);
    }

    handle_subblossom_shift(blossom, out_blossom);
    blossom->base = out_vertex;
    pathB.splice(pathB.end(), pathA);
    blossom->sub_blossoms = pathB; // pathB + pathA    
}

std::vector<BlossomMaximumMatching::Blossom*> 
BlossomMaximumMatching::blossoms_containing(NetworKit::node vertex, Blossom* until) {
    std::vector<Blossom*> res;
    Blossom* iter = trivial_blossom[vertex];

    while (iter != until) {
        res.push_back(iter);
        iter = iter->parent;
    }

    return res;
}

std::pair<
        std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::EdgeInfo>>, 
        std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::EdgeInfo>>>
BlossomMaximumMatching::split_subblossoms(
        std::list<std::pair<BlossomMaximumMatching::Blossom*, BlossomMaximumMatching::EdgeInfo>> sub_blossoms, 
        Blossom* blossom) {
    
    auto iter = sub_blossoms.begin();
    while (iter->first != blossom) iter ++;
    return { { sub_blossoms.begin(), std::next(iter) }, { std::next(iter), sub_blossoms.end() } };
}

// bool BlossomMaximumMatching::is_useful(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge) {
//     auto u_blossom = get_blossom(u);
//     auto v_blossom = get_blossom(v);
//     return !is_in_matching[edge] && edge_dual_variable(edge) == 0 && 
//             u_blossom->label == even && 
//             ((v_blossom->label == even && u_blossom != v_blossom) || v_blossom->label == free);
// }

bool BlossomMaximumMatching::is_exposed(BlossomMaximumMatching::Blossom* b) {
    return matched_vertex[b->base] == NetworKit::none;
}

void BlossomMaximumMatching::adjust_dual_variables() {
    auto delta1 = calc_delta1(); // min u_i : i - even vertex
    auto delta2 = calc_delta2(); // min pi_ij : i - even vertex, j - free vertex
    auto delta3 = calc_delta3(); // min pi_ij / 2 : i,j - even vertices in different blossoms
    auto delta4 = calc_delta4(); // min z_k : B_k - odd blossom 

    auto delta = std::min(delta1, std::min(delta2, std::min(delta3, delta4)));

    std::cerr << "delta_1 = " << delta1 << std::endl;
    std::cerr << "delta_2 = " << delta2 << std::endl;
    std::cerr << "delta_3 = " << delta3 << std::endl;
    std::cerr << "delta_4 = " << delta4 << std::endl;
    std::cerr << "Adjusting dual variables by " << delta << std::endl;
    adjust_by_delta(delta);

    if (delta == delta1) {
        // Optimal solution has been reached
        std::cerr << "Optimal solution reached\n";
        hasRun = true;
    } 
    else if (delta == delta3) find_delta3_useful_edges();
    else if (delta == delta2) find_delta2_useful_edges();
    else if (delta == delta4) {
        // Blossoms where minimum was achieved need to be expanded
        std::vector<Blossom*> to_expand;
        for (auto b : blossoms) {
            if (b->label == odd && b->z == 0.0 && !b->is_trivial()) {
                to_expand.push_back(b);
            }
        }
        for (auto b : to_expand) expand_odd_blossom(b);
    }
}

void BlossomMaximumMatching::print_backtrack(
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

void BlossomMaximumMatching::Blossom::delete_all_children() {
    for (auto [b, e] : sub_blossoms) {
        b->delete_all_children();
        delete b;
    }
}

BlossomMaximumMatching::Blossom::~Blossom() {
    // std::cerr << "Deleteing "; short_print(); std::cerr << std::endl;
    if (data != nullptr) delete data;
}

bool BlossomMaximumMatching::Blossom::is_trivial() {
    return sub_blossoms.size() == 0;
}

void BlossomMaximumMatching::Blossom::for_nodes(const std::function<void(NetworKit::node)>& handle) {
    if (is_trivial()) {
        handle(base);
    } else {
        for (auto [b, edge] : sub_blossoms)
            b->for_nodes(handle);
    }
}

void BlossomMaximumMatching::Blossom::check_consistency() {
    for (auto [b, e] : sub_blossoms) {
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

bool BlossomMaximumMatching::Blossom::contains(NetworKit::node v) {
    if (is_trivial()) return base == v;

    for (auto [b, e] : sub_blossoms) {
        if (b->contains(v)) return true;
    }
    return false;
}

void BlossomMaximumMatching::Blossom::print(int depth = 0) {
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
    for (auto [b, edge] : sub_blossoms) {
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
        for (auto [b, e] : sub_blossoms) {
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
        for (auto sb : sub_blossoms) sb.first->nodes_print();
    }
    std::cerr << ")";
}

BlossomMaximumMatching::EdgeInfo BlossomMaximumMatching::reverse(const EdgeInfo& info) {
    auto [u, v, id] = info;
    return { v, u, id };
}

} /* namespace Koala */
