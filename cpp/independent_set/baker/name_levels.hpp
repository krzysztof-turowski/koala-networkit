//
// Created by mikolajtwarog on 2021-04-30.
//

#ifndef TECHNIKA_BAKER_NAME_LEVELS_HPP
#define TECHNIKA_BAKER_NAME_LEVELS_HPP

int name_levels(PlanarEmbedding& embedding, std::vector<int>& outer_face, std::vector<int>& vertex_level,
                std::vector< std::vector<Edge> >& outer_edges) {
    for (int &v : vertex_level) {
        v = -1;
    }

    for (int v : outer_face) {
        vertex_level[v] = 1;
    }

    std::queue<Edge> next_level_edges;

    outer_edges.emplace_back();
    outer_edges.emplace_back();

    for (int i = 0; i < outer_face.size(); i++) {
        int v = outer_face[i];
        int w = outer_face[(i + 1) % outer_face.size()];
        outer_edges[1].push_back(Edge(v, w, nullptr));
        for (Edge e : embedding[v]) {
            if (vertex_level[e.m_source] == -1 || vertex_level[e.m_target] == -1) {
                next_level_edges.push(e);
            }

            if (e.m_source == w) {
                std::swap(e.m_source, e.m_target);
                int e_it = get_edge_it(e, w, embedding);
                std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
            }
        }
    }

    Edge next_level_edge;
    Edge current_edge;
    int starting_v;
    int current_v;
    int current_egde_it;

    int level = 1;
    std::vector<int> current_level;

    while (!next_level_edges.empty()) {
        current_level.clear();

        next_level_edge = next_level_edges.front();
        next_level_edges.pop();

        if (vertex_level[next_level_edge.m_source] > -1 && vertex_level[next_level_edge.m_target] > -1) {
            continue;
        }

        level = vertex_level[next_level_edge.m_source] == -1 ?
                vertex_level[next_level_edge.m_target] : vertex_level[next_level_edge.m_source];
        level++;

        starting_v = vertex_level[next_level_edge.m_source] == -1 ?
                     next_level_edge.m_source : next_level_edge.m_target;

        current_level.push_back(starting_v);

        current_egde_it = get_edge_it(next_level_edge, starting_v, embedding);

        vertex_level[starting_v] = level;
        bool res = false;
        for (int j = current_egde_it + 1; j < current_egde_it + embedding[starting_v].size(); j++) {
            Edge e_j = embedding[starting_v][j];
            if (vertex_level[e_j.m_source] == -1
                || vertex_level[e_j.m_target] == -1) {
                current_edge = e_j;
                res = true;
                break;
            }
        }
        if (!res) {
            continue;
        }

        vertex_level[starting_v] = -1;

        current_v = current_edge.m_source == starting_v ? current_edge.m_target : current_edge.m_source;
        current_egde_it = get_edge_it(current_edge, current_v, embedding);

        current_level.push_back(current_v);

        int next_to_starting_v = current_v;

        if (level >= outer_edges.size()) {
            outer_edges.emplace_back();
        }

        while (true) {
            outer_edges[level].push_back(current_edge);

            int temp_v = current_v;

            for (int j = current_egde_it + 1; j < current_egde_it + embedding[current_v].size(); j++) {
                Edge e_j = embedding[current_v][j];
                int neighbour = e_j.m_source == current_v ? e_j.m_target : e_j.m_source;
                if (vertex_level[neighbour] == -1) {
                    current_edge = e_j;
                    break;
                }
            }

            current_v = current_edge.m_source == current_v ? current_edge.m_target : current_edge.m_source;
            current_egde_it = get_edge_it(current_edge, current_v, embedding);

            current_level.push_back(current_v);

            if (temp_v == starting_v && current_v == next_to_starting_v) {
                break;
            }
        }

        current_level.pop_back();

        for (int i = 0; i < current_level.size(); i++) {
            int v = current_level[i];
            vertex_level[v] = level;
        }

        for (int i = 0; i < current_level.size(); i++) {
            int v = current_level[i];
            int w = current_level[(i+1) % current_level.size()];
            for (Edge& e : embedding[v]) {
                if (vertex_level[e.m_source] == -1
                    || vertex_level[e.m_target] == -1) {
                    next_level_edges.push(e);
                }

                if (e.m_source == w) {
                    std::swap(e.m_source, e.m_target);
                    int e_it = get_edge_it(e, w, embedding);
                    std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
                }
            }
        }
    }

    return level;
}

#endif //TECHNIKA_BAKER_NAME_LEVELS_HPP
