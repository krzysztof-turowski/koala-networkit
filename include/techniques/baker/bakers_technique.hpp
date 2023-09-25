#pragma once

#include <boost/graph/connected_components.hpp>

#include <techniques/baker/BakerKOuterplanar.hpp>
#include <techniques/baker/Bodlaender.hpp>

template <typename Problem, template<typename T> typename Algorithm>
int bakers_technique(Graph& g, NetworKit::Graph &G, PlanarEmbedding& embedding, std::vector<int>& outer_face, int k) {
    int res = Problem::getInitial();
    for (int r = 0; r <= Problem::getLimit(k); r++) {
        std::vector<int> vertex_level(num_vertices(g), -1);
        auto [max_level, _] = name_levels(embedding, outer_face, vertex_level);
        int num_of_graphs = Problem::getNumberOfGraphs(k, max_level);

        std::vector<std::vector<int>> levels(num_of_graphs);
        Problem::getLevels(vertex_level, r, k, levels);
        std::vector<std::vector<int>> level_vertices(num_vertices(g));
        for (int i = 0; i < levels.size(); i++) {
            for (int v : levels[i]) {
                level_vertices[v].push_back(i);
            }
        }
        std::vector<std::map<int, int>> global_to_local(num_of_graphs);
        for (int i = 0; i < levels.size(); i++) {
            for (int j = 0; j < levels[i].size(); j++) {
                global_to_local[i][levels[i][j]] = j;
            }
        }

        std::vector<PlanarEmbedding> level_embeddings(num_of_graphs);
        for (int i = 0; i < num_of_graphs; i++) {
            level_embeddings[i].resize(levels[i].size());
        }
        std::vector<Graph> level_graphs(num_of_graphs);
        std::vector<std::map<std::pair<int, int>, Edge>> added_edges(num_of_graphs);
        for (int v = 0; v < embedding.size(); v++) {
            if (level_vertices[v].empty()) {
                continue;
            }
            auto &edges = embedding[v];
            for (auto &e : edges) {
                if (level_vertices[e.m_source].empty() || level_vertices[e.m_target].empty()) {
                    continue;
                }

                std::vector<int> intersection;
                std::set_intersection(level_vertices[e.m_source].begin(), level_vertices[e.m_source].end(),
                                      level_vertices[e.m_target].begin(), level_vertices[e.m_target].end(),
                                      std::back_inserter(intersection));

                for (int level : intersection) {
//                int level = intersection.front();
                    int local_source = global_to_local[level][e.m_source];
                    int local_target = global_to_local[level][e.m_target];
                    std::pair<int, int> e_pair(std::min(e.m_source, e.m_target), std::max(e.m_source, e.m_target));
                    auto local_e = added_edges[level].find(e_pair);
                    if (local_e == added_edges[level].end()) {
                        Edge new_e = add_edge(local_source, local_target, level_graphs[level]).first;
                        level_embeddings[level][global_to_local[level][v]].push_back(new_e);
                        added_edges[level][e_pair] = new_e;
                    } else {
                        level_embeddings[level][global_to_local[level][v]].push_back(local_e->second);
                    }
                }
            }
        }

        int temp_res = 0;

        for (int i = 0; i < level_graphs.size(); i++) {
            Graph &level_g = level_graphs[i];
            if (num_vertices(level_g) == 0) {
                res++;
                continue;
            }

            std::vector<int> v_comp(num_vertices(level_g));
            int num = connected_components(level_g, &v_comp[0]);

            std::vector<std::vector<int> > components(num);

            for (int v = 0; v < num_vertices(level_g); v++) {
                components[v_comp[v]].push_back(v);
            }

            global_to_local.clear();
            global_to_local.resize(num);

            for (int c = 0; c < components.size(); c++) {
                for (int j = 0; j < components[c].size(); j++) {
                    global_to_local[c][components[c][j]] = j;
                }
            }

            std::vector<PlanarEmbedding> conn_embeddings(num);
            for (int j = 0; j < num; j++) {
                conn_embeddings[j].resize(components[j].size());
            }
            std::vector<Graph> conn_graphs(num);
            added_edges.clear();
            added_edges.resize(num);

            for (int v = 0; v < level_embeddings[i].size(); v++) {
                for (auto &e : level_embeddings[i][v]) {
                    if (v_comp[e.m_source] == v_comp[e.m_target]) {
                        int comp = v_comp[e.m_source];
                        int local_source = global_to_local[comp][e.m_source];
                        int local_target = global_to_local[comp][e.m_target];
                        std::pair<int, int> e_pair(std::min(e.m_source, e.m_target), std::max(e.m_source, e.m_target));
                        auto local_e = added_edges[comp].find(e_pair);
                        if (local_e == added_edges[comp].end()) {
                            Edge new_e = add_edge(local_source, local_target, conn_graphs[comp]).first;
                            conn_embeddings[comp][global_to_local[comp][v]].push_back(new_e);
                            added_edges[comp][e_pair] = new_e;
                        } else {
                            conn_embeddings[comp][global_to_local[comp][v]].push_back(local_e->second);
                        }
                    }
                }
            }

            for (int c = 0; c < num; c++) {
                Graph &graph = conn_graphs[c];
                PlanarEmbedding &emb = conn_embeddings[c];

                if (num_vertices(graph) <= 1) {
                    temp_res++;
                } else {
                    std::map<std::pair<int, int>, std::vector<int>> faces;
                    std::vector<std::vector<int>> vertices_in_face;
                    face_getter visitor(faces, vertices_in_face);
                    PlanarEmbedding2 embedding2(convert(embedding));
                    level_face_traversal(embedding2, visitor);

                    int min_level = INT16_MAX;
                    for (int v : components[c]) {
                        min_level = std::min(min_level, vertex_level[levels[i][v]]);
                    }

                    std::set<int> conn_outer_face;
                    for (int v : components[c]) {
                        if (vertex_level[levels[i][v]] == min_level) {
                            conn_outer_face.insert(global_to_local[c][v]);
                        }
                    }

                    std::vector<int> out_face;
                    for (auto &face : vertices_in_face) {
                        std::set<int> temp_face(face.begin(), face.end());
                        if (conn_outer_face == temp_face) {
                            out_face = std::vector<int>(face.begin(), face.end());
                            break;
                        }
                    }
                    temp_res += Algorithm<Problem>(graph, G, embedding, out_face);
                }
            }
        }
        res = Problem::getBest(res, temp_res);
    }
    return res;
}
