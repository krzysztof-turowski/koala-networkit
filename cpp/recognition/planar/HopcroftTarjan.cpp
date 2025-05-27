#include <iostream>
#include <list>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>
#include <recognition/planar/PlanarGraphRecognition.hpp>

namespace Koala {

    class Block {
    private:
        std::list<int> Latt, Ratt;
        std::list<NetworKit::WeightedEdge> Lseg, Rseg;

    public:
        Block(NetworKit::WeightedEdge e, std::list<int> &A) {
            Lseg.push_back(e);
            Latt.splice(Latt.end(), A);
        }

        ~Block() {}

        void flip() {
            std::list<int> ha;
            std::list<NetworKit::WeightedEdge> he;
            ha.splice(ha.end(), Ratt);
            Ratt.splice(Ratt.end(), Latt);
            Latt.splice(Latt.end(), ha);
            he.splice(he.end(), Rseg);
            Rseg.splice(Rseg.end(), Lseg);
            Lseg.splice(Lseg.end(), he);
        }

        bool left_interlace(std::stack<Block> &S) {
            return (!S.empty() && !S.top().Latt.empty() &&
                    Latt.back() < S.top().Latt.front());
        }

        bool right_interlace(std::stack<Block> &S) {
            return (!S.empty() && !S.top().Ratt.empty() &&
                    Latt.back() < S.top().Ratt.front());
        }

        void combine(Block &b) {
            Latt.splice(Latt.end(), b.Latt);
            Ratt.splice(Ratt.end(), b.Ratt);
            Lseg.splice(Lseg.end(), b.Lseg);
            Rseg.splice(Rseg.end(), b.Rseg);
        }

        bool clean(int dfsnum_w, std::vector<int> &alpha,
                   std::vector<NetworKit::node> &dfsnum) {
            while (!Latt.empty() && Latt.front() == dfsnum_w) {
                Latt.pop_front();
            }
            while (!Ratt.empty() && Ratt.front() == dfsnum_w) {
                Ratt.pop_front();
            }
            if (!Latt.empty() || !Ratt.empty()) {
                return false;
            }
            for (auto e: Lseg) {
                alpha[e.weight] = LEFT;
            }
            for (auto e: Rseg) {
                alpha[e.weight] = RIGHT;
            }
            return true;
        }

        void add_to_Att(std::list<int> &Att, int dfsnum_w0,
                        std::vector<int> &alpha) {
            if (!Ratt.empty() && this->Ratt.front() > dfsnum_w0) {
                flip();
            }
            Att.splice(Att.end(), Latt);
            Att.splice(Att.end(), Ratt);
            for (auto e: Lseg) {
                alpha[e.weight] = LEFT;
            }
            for (auto e: Rseg) {
                alpha[e.weight] = RIGHT;
            }
        }

        friend class HopcroftTarjan;
    };

    HopcroftTarjan::HopcroftTarjan(NetworKit::Graph &graph, bool embedding)
            : PlanarGraphRecognition(graph, embedding) {}

    void HopcroftTarjan::dfs(NetworKit::Graph &g, NetworKit::node v,
                             std::vector<bool> &reached) {
        reached[v] = true;
        for (auto u: g.neighborRange(v)) {
            if (!reached[u]) {
                dfs(g, u, reached);
            }
        }
    }

    void HopcroftTarjan::dfs_in_make_biconnected_graph(
            NetworKit::Graph &g, NetworKit::node v, int &dfs_count,
            std::vector<bool> &reached, std::vector<NetworKit::node> &dfsnum,
            std::vector<NetworKit::node> &lowpt, std::vector<NetworKit::node> &parent) {
        NetworKit::node w = NetworKit::none;
        dfsnum[v] = dfs_count++;
        lowpt[v] = dfsnum[v];
        reached[v] = true;
        std::vector<int> neighbors;
        for (auto u: g.neighborRange(v)) {
            neighbors.push_back(u);
        }
        for (auto u: neighbors) {
            if (!reached[u]) {
                if (w == NetworKit::none) {
                    w = u;
                }
                parent[u] = v;
                dfs_in_make_biconnected_graph(g, u, dfs_count, reached, dfsnum,
                                              lowpt, parent);
                if (lowpt[u] == dfsnum[v]) {
                    if (u == w && parent[v] != NetworKit::none) {
                        g.addEdge(u, parent[v]);
                        g.addEdge(parent[v], u);
                    }
                    if (u != w) {
                        g.addEdge(u, w);
                        g.addEdge(w, u);
                    }
                }
                lowpt[v] = std::min(lowpt[v], lowpt[u]);
            } else {
                lowpt[v] = std::min(lowpt[v], dfsnum[u]);
            }
        }
    }

    void HopcroftTarjan::make_biconnected_graph(NetworKit::Graph &g) {
        // make it connected
        NetworKit::node u = 0;
        std::vector<bool> reached(N);
        for (auto v: g.nodeRange()) {
            if (!reached[v]) {
                dfs(g, v, reached);
                if (u != v) {
                    g.addEdge(u, v);
                    g.addEdge(v, u);
                }
            }
        }
        // make it biconnected
        fill(reached.begin(), reached.end(), false);
        std::vector<NetworKit::node> lowpt(N);
        std::vector<NetworKit::node> parent(N, NetworKit::none);
        std::vector<NetworKit::node> dfsnum(N);
        int dfs_count = 0;

        dfs_in_make_biconnected_graph(g, u, dfs_count, reached, dfsnum, lowpt,
                                      parent);
    }

    void HopcroftTarjan::dfs_in_reorder(NetworKit::Graph &g, NetworKit::node v,
                                        std::vector<NetworKit::Edge> &Add,
                                        int &dfs_count, std::vector<bool> &reached,
                                        std::vector<NetworKit::node> &dfsnum,
                                        std::vector<NetworKit::node> &lowpt,
                                        std::vector<NetworKit::node> &lowpt2,
                                        std::vector<NetworKit::node> &parent) {
        dfsnum[v] = dfs_count++;
        lowpt2[v] = lowpt[v] = dfsnum[v];
        reached[v] = true;
        for (auto w: g.neighborRange(v)) {
            if (!reached[w]) {
                Add.push_back(NetworKit::Edge(v, w));
                parent[w] = v;
                dfs_in_reorder(g, w, Add, dfs_count, reached, dfsnum, lowpt, lowpt2,
                               parent);
                lowpt[v] = std::min(lowpt[v], lowpt[w]);
            } else {
                lowpt[v] = std::min(lowpt[v], dfsnum[w]);
                if (dfsnum[w] < dfsnum[v] && parent[v] != w) {
                    Add.push_back({v, w});
                }
            }
        }
        for (auto w: g.neighborRange(v)) {
            if (parent[w] == v) {
                if (lowpt[w] != lowpt[v]) {
                    lowpt2[v] = std::min(lowpt2[v], lowpt[w]);
                }
                lowpt2[v] = std::min(lowpt2[v], lowpt2[w]);
            } else {
                if (lowpt[v] != dfsnum[w]) {
                    lowpt2[v] = std::min(lowpt2[v], dfsnum[w]);
                }
            }
        }
    }

    void HopcroftTarjan::reorder(NetworKit::Graph &g,
                                 std::vector<NetworKit::node> &parent,
                                 std::vector<NetworKit::node> &dfsnum,
                                 std::vector<NetworKit::Edge> &order) {
        NetworKit::node v = 0;
        std::vector<bool> reached(N);
        int dfs_count = 0;
        std::vector<NetworKit::Edge> Add;
        std::vector<NetworKit::node> lowpt(N), lowpt2(N);
        dfs_in_reorder(g, v, Add, dfs_count, reached, dfsnum, lowpt, lowpt2,
                       parent);
        std::vector<NetworKit::Edge> buckets[N * 2 + 5];
        for (auto[u, v]: Add) {
            int cost = ((dfsnum[v] < dfsnum[u])
                        ? 2 * dfsnum[v]
                        : ((lowpt2[v] >= dfsnum[u]) ? 2 * lowpt[v]
                                                    : 2 * lowpt[v] + 1));
            buckets[cost].push_back(NetworKit::Edge(u, v));
        }
        g = NetworKit::Graph(N, true, true);
        order = std::vector<NetworKit::Edge>(Add.size());
        int id = 0;
        for (int i = 0; i <= 2 * N; i++) {
            for (auto[u, v]: buckets[i]) {
                g.addEdge(u, v, id);
                order[id] = NetworKit::Edge(u, v);
                id++;
            }
        }
    }

    bool HopcroftTarjan::strong_planar(NetworKit::Graph &g, NetworKit::Edge e,
                                       std::list<int> &Att, std::vector<int> &alpha,
                                       std::vector<NetworKit::node> &dfsnum,
                                       std::vector<NetworKit::node> &parent) {
        NetworKit::node y = e.v, t = *g.neighborRange(y).begin(), x = e.u, wk = y;
        while (dfsnum[t] > dfsnum[wk]) {
            wk = t;
            t = *g.neighborRange(wk).begin();
        }
        NetworKit::node w0 = t;
        NetworKit::node w = wk;
        std::stack<Block> S;
        while (w != x) {
            int count = 0;
            for (auto[v, weight]: g.weightNeighborRange(w)) {
                count++;
                if (count != 1) {
                    std::list<int> A;
                    if (dfsnum[w] < dfsnum[v]) {
                        if (!strong_planar(g, NetworKit::Edge(w, v), A, alpha,
                                           dfsnum, parent)) {
                            while (!S.empty()) {
                                S.pop();
                            }
                            return false;
                        }
                    } else {
                        A.push_back(dfsnum[v]);
                    }

                    Block b(NetworKit::WeightedEdge(w, v, weight), A);
                    while (true) {
                        if (b.left_interlace(S)) {
                            S.top().flip();
                        }
                        if (b.left_interlace(S)) {
                            while (!S.empty()) {
                                S.pop();
                            }
                            return false;
                        }
                        if (b.right_interlace(S)) {
                            b.combine(S.top());
                            S.pop();
                        } else {
                            break;
                        }
                    }
                    S.push(b);
                }
            }
            while (!S.empty() &&
                   S.top().clean(dfsnum[parent[w]], alpha, dfsnum)) {
                S.pop();
            }
            w = parent[w];
        }
        Att.clear();
        while (!S.empty()) {
            Block b = S.top();
            S.pop();
            if (!(b.Latt.empty()) && !(b.Ratt.empty()) &&
                (b.Latt.front() > dfsnum[w0]) &&
                (b.Ratt.front() > dfsnum[w0])) {
                while (!S.empty()) {
                    S.pop();
                }
                return false;
            }
            b.add_to_Att(Att, dfsnum[w0], alpha);
        }
        if (w0 != x) {
            Att.push_back(dfsnum[w0]);
        }

        return true;
    }

    PlanarGraphRecognition::State HopcroftTarjan::planar() {
        N = graph.numberOfNodes();
        if (N <= 3) {
            return PlanarGraphRecognition::State::PLANAR;
        }
        if (graph.numberOfEdges() > 6 * N - 12) {
            return PlanarGraphRecognition::State::NOT_PLANAR;
        }

        // make G copy of Graph and make it biconnected

        NetworKit::Graph g = NetworKit::Graph(N, false, true);
        for (auto i: graph.edgeRange()) {
            g.addEdge(i.u, i.v);
            g.addEdge(i.v, i.u);
        }
        make_biconnected_graph(g);
        if (g.numberOfEdges() > 6 * N - 12) {
            return PlanarGraphRecognition::State::NOT_PLANAR;
        }

        // test planarity
        std::vector<NetworKit::node> parent(N, NetworKit::none);
        std::vector<NetworKit::node> dfsnum(N);
        std::vector<NetworKit::Edge> order;
        reorder(g, parent, dfsnum, order);

        std::vector<int> alpha(g.numberOfEdges());
        std::list<int> Att;
        NetworKit::Edge e(0, *g.neighborRange(0).begin());
        alpha[0] = LEFT;

        if (!strong_planar(g, e, Att, alpha, dfsnum, parent)) {
            return PlanarGraphRecognition::State::NOT_PLANAR;
        }

        if (embedding) {
            // TODO embedding
        }
        return PlanarGraphRecognition::State::PLANAR;
    }

    void HopcroftTarjan::run() {
        hasRun = true;
        this->is_planar = this->planar();
    }

}  // namespace Koala
