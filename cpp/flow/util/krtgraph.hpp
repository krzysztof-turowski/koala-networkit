#pragma once

#include <iostream>
#include <utility>
#include <vector>
#include <bits/stdc++.h>

#define PII std::pair<int, int>

struct KRTEdge {
    int dest, cost, flow, backLink;

    explicit KRTEdge(int d, int c = 0, int f = 0, int bl = -1) :
            dest(d), cost(c), flow(f), backLink(bl) {}
};

class KRTGraph {
private:
    int n{}, m{}, s{}, t{};
    std::vector<std::vector<KRTEdge>> edges;

    KRTGraph(int, int, int, int, std::map<std::pair<int, int>, int>);

public:
    KRTGraph();

    explicit KRTGraph(std::ifstream &);

    KRTGraph(int, int, int, int, std::vector<std::vector<KRTEdge>>);

    int getN() const;

    int getM() const;

    int getS() const;

    int getT() const;

    const std::vector<std::vector<KRTEdge>> &getEdges() const;

    static KRTGraph read_from_file(const char *);

    static KRTGraph read_from_stdin();
};

KRTGraph::KRTGraph() : n(0), m(0), s(0), t(0) {}

KRTGraph::KRTGraph(std::ifstream &in) {
    in.read(reinterpret_cast<char *>(&n), sizeof(n));
    in.read(reinterpret_cast<char *>(&m), sizeof(m));
    in.read(reinterpret_cast<char *>(&s), sizeof(s));
    in.read(reinterpret_cast<char *>(&t), sizeof(t));

    edges = std::vector<std::vector<KRTEdge>>(n + 1);

    std::vector<int> tmpIn(3 * m);
    in.read(reinterpret_cast<char *>(tmpIn.data()), 3 * m * sizeof(int));

    std::map<PII, int> edgesList;
    for (std::size_t i = 0; i < tmpIn.size(); i += 3) {
        edgesList[std::make_pair(tmpIn[i], tmpIn[i + 1])] += tmpIn[i + 2];
    }

    *this = KRTGraph(n, m, s, t, edgesList);
}

KRTGraph::KRTGraph(int n, int m, int s, int t, std::vector<std::vector<KRTEdge>> edges)
        : n(n), m(m), s(s), t(t), edges(std::move(edges)) {}

KRTGraph KRTGraph::read_from_file(const char *filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return KRTGraph();
    }
    return KRTGraph(in);
}

KRTGraph KRTGraph::read_from_stdin() {
    int n, m, s, t;
    std::cin >> n >> m >> s >> t;

    std::map<PII, int> edgesList;
    for (int i = 0; i < m; ++i) {
        int a, b, c;
        std::cin >> a >> b >> c;
//        std::cerr << a << " " << b << " " << c << std::endl;
        edgesList[std::make_pair(a, b)] += c;
    }
//    std::cerr << edgesList.size() << std::endl;
    return KRTGraph(n, m, s, t, edgesList);
}

int KRTGraph::getN() const {
    return n;
}

int KRTGraph::getM() const {
    return m;
}

int KRTGraph::getS() const {
    return s;
}

int KRTGraph::getT() const {
    return t;
}

const std::vector<std::vector<KRTEdge>> &KRTGraph::getEdges() const {
    return edges;
}

KRTGraph::KRTGraph(int n, int m, int s, int t, std::map<std::pair<int, int>, int> edgesList) :
        n(n), m(m), s(s), t(t) {
    edges = std::vector<std::vector<KRTEdge>>(n + 1);
    for (auto &x : edgesList) {
        if (x.second < 0) continue;

        PII info = x.first;
        int val = x.second;

        edges[info.first].push_back(KRTEdge(info.second, val, 0, edges[info.second].size()));

        // TODO: Uncomment and handle properly in PR
//        auto rev = edgesList.find(std::make_pair(info.second, info.first));
//        if (rev != edgesList.end()) {
//            PII revInfo = rev->first;
//            int revVal = rev->second;
//            edges[revInfo.first].push_back(edge(revInfo.second, revVal, 0, edges[revInfo.second].size() - 1));
//            rev->second = -1;
//        } else {
//            edges[info.second].push_back(edge(info.first, 0, 0, edges[info.first].size() - 1));
//        }
    }
}
