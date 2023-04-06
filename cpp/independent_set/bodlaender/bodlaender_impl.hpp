//
// Created by mikolajtwarog on 2021-05-24.
//

#ifndef TECHNIKA_BAKER_BODLAENDER_IMPL_HPP
#define TECHNIKA_BAKER_BODLAENDER_IMPL_HPP

#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include "../baker/create_tree_decomposition.hpp"

using namespace boost;
namespace mp = multiprecision;

enum Problem {vc, is, ds_ecc, ds_lcc};

class bodlaender_impl {

    struct tr_node{
        dynamic_bitset<> f;
        std::vector<int> f2;
        std::vector<int> r_values;
        int x;

        tr_node(int n): f(n, 0), f2(n, -1) {}
        tr_node(int n, int f_num): f(n, f_num) {}
    };

    tree_decomposition tr;
    Graph graph;
    PlanarEmbedding embedding;
    std::vector< std::vector<tr_node> > tables;
    int n;
    Problem problem;
    bool minimum;


    std::vector< std::set<int> > incidence;

    void calculate_table(int v, int p) {
        std::vector<int>& children = tr.tree[v];

        for (int child : children) {
            if (child != p) {
                calculate_table(child, v);
            }
        }


        std::vector<tr_node> temp;
        std::set<int> vertices;
        vertices = std::set<int>(tr.nodes[v].begin(), tr.nodes[v].end());
        Graph &sub_g = graph.create_subgraph(vertices.begin(), vertices.end());

        if (problem == ds_ecc) {
//            for (int g_v : vertices) {
//                int l_g_v = sub_g.global_to_local(g_v);
//                typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
//                for (boost::tie(ei, ei_end) = out_edges(l_g_v, sub_g); ei != ei_end; ++ei) {
//                    int nei = sub_g.local_to_global(l_g_v == ei->m_source ? ei->m_target : ei->m_source);
//                    incidence[g_v].insert(nei);
//                }
//            }

            mp::cpp_int num = mp::pow(mp::cpp_int(n), vertices.size());

            for (mp::cpp_int x = 0; x < num; x++) {
                temp.emplace_back(n);
                auto &f = temp.back().f2;
                mp::cpp_int x_temp = x;
                int res = 0;
                for (int g_v : vertices) {
                    if (int(x_temp % n) == g_v) {
                        res++;
                    }
                    f[g_v] = int(x_temp % n);
                    x_temp /= n;
                }

                int r_1 = 1, r_2 = 1, r_3 = 1;

                typename graph_traits<Graph>::edge_iterator ei, ei_end;
                for (boost::tie(ei, ei_end) = edges(sub_g); ei != ei_end; ++ei) {
                    int v_1 = sub_g.local_to_global(ei->m_source), v_2 = sub_g.local_to_global(ei->m_target);
                    r_1 &= f[v_1] != v_2 || f[v_2] == v_2;
                    r_2 &= f[v_2] != v_1 || f[v_1] == v_1;
                }

                for (int g_v : vertices) {
                    r_3 &= incidence[g_v].find(f[g_v]) != incidence[g_v].end() || f[g_v] == g_v;
                }

//                temp.back().r_values.push_back(r_1);
//                temp.back().r_values.push_back(r_2);
//                temp.back().r_values.push_back(r_3);
                temp.back().r_values.push_back(r_1 & r_2 & r_3);
                temp.back().r_values.push_back(res);
            }
        } else {
            int num = 1 << vertices.size();
            for (int x = 0; x < num; x++) {
                temp.emplace_back(n);
                auto &f = temp.back().f;
                int x_temp = x;
                for (int g_v : vertices) {
                    if (x_temp & 1) {
                        f[g_v] = 1;
                    }
                    x_temp >>= 1;
                }

                if (problem == ds_lcc) {
                    int r_val = 1;
                    for (int g_v : vertices) {
                        int r = f[g_v];
                        for (Edge e : embedding[g_v]) {
                            int nei = g_v == e.m_source ? e.m_target : e.m_source;
                            if (vertices.find(nei) != vertices.end()) {
                                r += f[nei];
                            } else {
                                r = 1;
                                break;
                            }
                        }
                        r_val &= (r > 0);
                    }
                    temp.back().r_values.push_back(r_val);
                } else {
                    int r_val = 1;
                    for (int g_v : vertices) {
                        int l_g_v = sub_g.global_to_local(g_v);
                        typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
                        for (boost::tie(ei, ei_end) = out_edges(l_g_v, sub_g); ei != ei_end; ++ei) {
                            int nei = sub_g.local_to_global(l_g_v == ei->m_source ? ei->m_target : ei->m_source);
                            if (vertices.find(nei) != vertices.end()) {
                                if (minimum) {
                                    r_val &= (f[g_v] + f[nei] > 0);
                                } else {
                                    r_val &= (f[g_v] + f[nei] < 2);
                                }
                            }
                        }
                    }
                    temp.back().r_values.push_back(r_val);
                }

                temp.back().r_values.push_back(f.count());
                temp.back().x = x;
            }
        }

        for (int child : children) {
            if (child != p) {
                temp = calculate_temp(v, child, temp);
            }
        }
        tables[v] = temp;
    }

    std::vector<tr_node> calculate_temp (int v, int child, std::vector<tr_node>& temp) {
        std::vector<tr_node> new_temp;
        std::vector<int> intersection;
        std::set_intersection(tr.nodes[v].begin(), tr.nodes[v].end(), tr.nodes[child].begin(),
                              tr.nodes[child].end(), std::back_inserter(intersection));

        for (auto& f1 : temp) {
            for (auto& f2 : tables[child]) {
                bool good = true;
                for (int ver : intersection) {
                    if (f1.f[ver] != f2.f[ver] || f1.f2[ver] != f2.f2[ver]) {
                        good = false;
                        break;
                    }
                }

                if (!good) {
                    continue;
                }

                std::vector<int> s_values;
                for (int i = 0; i < f1.r_values.size() - 1; i++) {
                    s_values.push_back(f1.r_values[i] & f2.r_values[i]);
                }

                if (problem == ds_ecc) {
                    int r_m = 0;
                    for (int i = 0; i < n ; i++) {
                        if (f1.f2[i] == i && f2.f2[i] == i) {
                            r_m++;
                        }
                    }
                    s_values.push_back(f1.r_values.back() + f2.r_values.back() - r_m);
                } else {
                    auto f3 = f1.f & f2.f;
                    s_values.push_back(f1.r_values.back() + f2.r_values.back() - f3.count());
                }

                int f_it = -1;
                for (int i = 0; i < new_temp.size(); i++) {
                    good = true;
                    for (int r = 0; r < s_values.size() - 1; r++) {
                        if (s_values[r] != new_temp[i].r_values[r]) {
                            good = false;
                            break;
                        }
                    }

                    if (good && f1.f == new_temp[i].f && (problem != ds_ecc || f1.f2 == new_temp[i].f2)) {
                        f_it = i;
                        break;
                    }
                }

                if (f_it == -1) {
                    new_temp.emplace_back(n);
                    new_temp.back().f = f1.f;
                    if (problem == ds_ecc) {
                        new_temp.back().f2 = f1.f2;
                    }
                    new_temp.back().r_values = s_values;
                    new_temp.back().x = f1.x | f2.x;
                } else {
                    if ((minimum && new_temp[f_it].r_values.back() > s_values.back())
                    || (!minimum && new_temp[f_it].r_values.back() < s_values.back())) {
                        new_temp[f_it].r_values.back() = s_values.back();
                        new_temp[f_it].f = f1.f;
                        if (problem == ds_ecc) {
                            new_temp[f_it].f2 = f1.f2;
                        }
                        new_temp[f_it].x = f1.x | f2.x;
                    }
                }
            }
        }

        return new_temp;
    }

public:
    bodlaender_impl(Graph g, PlanarEmbedding emb, std::vector<int> outer_face, Problem pr):
    graph(g), embedding(emb), n(num_vertices(g)), problem(pr) {
        minimum = pr == vc || pr == ds_ecc || pr == ds_lcc;
        get_tree_decomposition(g, emb, outer_face, tr);
        tables.resize(tr.tree.size());

        if (problem == ds_ecc) {
            for (int v = 0; v < embedding.size(); v++) {
                incidence.emplace_back();

                for (auto& e : embedding[v]) {
                    int w = v == e.m_source ? e.m_target : e.m_source;
                    incidence.back().insert(w);
                }
            }
        }

        if (problem == ds_lcc) {
            for (auto &node : tr.nodes) {
                std::vector<int> to_add;
                for (int v : node) {
                    for (auto &e : embedding[v]) {
                        int w = v == e.m_source ? e.m_target : e.m_source;
                        to_add.push_back(w);
                    }
                }

                for (int v : to_add) {
                    node.insert(v);
                }
            }
        }

        calculate_table(0, -1);
    }

    int get_result() {
        int result;
        if (minimum) {
            result = INT16_MAX;
        } else {
            result = -1;
        }

        for (auto& node : tables[0]) {
            int good = true;
            for (int r = 0; r < node.r_values.size() - 1; r++) {
                if (node.r_values[r] == 0) {
                    good = false;
                    break;
                }
            }

            if (good) {
                if (minimum) {
                    result = std::min(result, node.r_values.back());
                } else {
                    result = std::max(result, node.r_values.back());
                }
            }
        }

        return result;
    }
};

int bodlaender_vertex_cover(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender_impl bd(g, emb, outer_face, vc);
    return bd.get_result();
}

int bodlaender_independent_set(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender_impl bd(g, emb, outer_face, is);
    return bd.get_result();
}

int bodlaender_dominating_set_lcc(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender_impl bd(g, emb, outer_face, ds_lcc);
    return bd.get_result();
}

int bodlaender_dominating_set_ecc(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender_impl bd(g, emb, outer_face, ds_ecc);
    return bd.get_result();
}

int bodlaender(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face, Problem p) {
    bodlaender_impl bd(g, emb, outer_face, p);
    return bd.get_result();
}


#endif //TECHNIKA_BAKER_BODLAENDER_IMPL_HPP
