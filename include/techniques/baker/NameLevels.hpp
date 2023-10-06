/*
 * NameLevels.hpp
 *
 *  Created on: 30.04.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

void compare(const PlanarEmbedding &emb, const PlanarEmbedding2 &emb2) {
    assert(emb.size() == emb2.size());
    for (int i = 0; i < emb.size(); i++) {
        assert(emb[i].size() == emb2[i].size());
        for (int j = 0; j < emb[i].size(); j++) {
          assert(emb[i][j].m_source == emb2[i][j].u);
          assert(emb[i][j].m_target == emb2[i][j].v);
        }
    }
}

PlanarEmbedding2 convert(const PlanarEmbedding &emb) {
    PlanarEmbedding2 emb2(emb.size());
    for (int i = 0; i < emb.size(); i++) {
        emb2[i].reserve(emb[i].size());
        for (int j = 0; j < emb[i].size(); j++) {
            emb2[i].emplace_back(emb[i][j].m_source, emb[i][j].m_target);
        }
    }
    compare(emb, emb2);
    return emb2;
}

PlanarEmbedding convert(const PlanarEmbedding2 &emb2) {
    PlanarEmbedding emb(emb2.size());
    for (int i = 0; i < emb2.size(); i++) {
        emb[i].reserve(emb2[i].size());
        for (int j = 0; j < emb2[i].size(); j++) {
            emb[i].emplace_back(emb2[i][j].u, emb2[i][j].v, nullptr);
        }
    }
    compare(emb, emb2);
    return emb;
}

auto name_levels(
        const PlanarEmbedding2& embedding, std::vector<int>& outer_face, std::vector<int>& vertex_level) {
    for (int &v : vertex_level) {
        v = -1;
    }
    for (int v : outer_face) {
        vertex_level[v] = 1;
    }
    std::vector<std::vector<NetworKit::Edge>> outer_edges(2);
    std::queue<NetworKit::Edge> next_level_edges;
    for (int i = 0; i < outer_face.size(); i++) {
        int v = outer_face[i], w = outer_face[(i + 1) % outer_face.size()];
        outer_edges[1].emplace_back(v, w);
        for (const auto &e : embedding[v]) {
            if (vertex_level[e.u] == -1 || vertex_level[e.v] == -1) {
                next_level_edges.emplace(e.u, e.v);
            }
        }
    }
    int level = 1;
    while (!next_level_edges.empty()) {
        std::vector<int> current_level;
        auto [u, v] = next_level_edges.front();
        next_level_edges.pop();
        if (vertex_level[u] > -1 && vertex_level[v] > -1) {
            continue;
        }
        level = std::max(vertex_level[v], vertex_level[u]) + 1;
        int starting_v = vertex_level[u] == -1 ? u : v;
        current_level.push_back(starting_v);
        int current_edge_it = get_edge_it(starting_v, starting_v == u ? v : u, embedding);
        vertex_level[starting_v] = level;
        std::optional<NetworKit::Edge> current_edge = std::nullopt;
        for (int j = current_edge_it + 1; j < current_edge_it + embedding[starting_v].size(); j++) {
            auto e_j = NetworKit::Edge(embedding[starting_v][j].u, embedding[starting_v][j].v);
            if (vertex_level[e_j.u] == -1 || vertex_level[e_j.v] == -1) {
                current_edge = std::make_optional(e_j);
                break;
            }
        }
        if (!current_edge.has_value()) {
            continue;
        }
        vertex_level[starting_v] = -1;
        int current_v = current_edge->u == starting_v ? current_edge->v : current_edge->u;
        current_edge_it = get_edge_it(*current_edge, current_v, embedding);
        current_level.push_back(current_v);
        int next_to_starting_v = current_v;
        if (level >= outer_edges.size()) {
            outer_edges.emplace_back();
        }
        while (true) {
            outer_edges[level].emplace_back(current_edge->u, current_edge->v);
            int temp_v = current_v;
            for (int j = current_edge_it + 1; j < current_edge_it + embedding[current_v].size(); j++) {
                auto e_j = NetworKit::Edge(embedding[current_v][j].u, embedding[current_v][j].v);
                int neighbour = e_j.u == current_v ? e_j.v : e_j.u;
                if (vertex_level[neighbour] == -1) {
                    current_edge = std::make_optional(e_j);
                    break;
                }
            }
            current_v = current_edge->u == current_v ? current_edge->v : current_edge->u;
            current_edge_it = get_edge_it(*current_edge, current_v, embedding);
            current_level.push_back(current_v);
            if (temp_v == starting_v && current_v == next_to_starting_v) {
                break;
            }
        }
        current_level.pop_back();
        for (int v : current_level) {
            vertex_level[v] = level;
        }
        for (int i = 0; i < current_level.size(); i++) {
            int v = current_level[i], w = current_level[(i + 1) % current_level.size()];
            for (const auto &e : embedding[v]) {
                if (vertex_level[e.u] == -1 || vertex_level[e.v] == -1) {
                    next_level_edges.emplace(e.u, e.v);
                }
            }
        }
    }
    return std::make_tuple(level, outer_edges);
}
