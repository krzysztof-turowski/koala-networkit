#pragma once

#include <map>
#include <set>
#include <list>
#include <optional>
#include <functional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <matching/PriorityQueues.hpp>

namespace Koala {

/**
 * @ingroup matching
 * The base class for the maximum matching algorithms.
 *
 */
class MaximumMatching : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the maximum matching algorithm.
     *
     * @param graph The input graph.
     */
    MaximumMatching(NetworKit::Graph &graph);

    /**
     * Return the matching found by the algorithm.
     *
     * @return a matching between nodes.
     */
    const std::map<NetworKit::node, NetworKit::node>& getMatching() const;

protected:
    NetworKit::Graph graph;
    std::map<NetworKit::node, NetworKit::node> matching;
};


/**
 * 
*/
class BlossomMaximumMatching : public MaximumMatching {
public:
    BlossomMaximumMatching(NetworKit::Graph &graph);

    /**
     * Execute the Edmonds maximum matching algorithm.
     */
    void run();

protected:
    virtual ~BlossomMaximumMatching();

    struct EdgeInfo { 
        NetworKit::node u, v; 
        NetworKit::edgeid id; 
    };

    static constexpr EdgeInfo no_edge {NetworKit::none, NetworKit::none,NetworKit::none};

    enum BlossomLabel { odd, even, free };

    class Blossom {
    public:
        class BlossomData {};
        Blossom* parent;
        NetworKit::node base;
        std::list<std::pair<Blossom*, EdgeInfo>> sub_blossoms;
        BlossomLabel label;
        EdgeInfo backtrack_edge;
        bool visited;
        NetworKit::edgeweight z;
        BlossomData* data;

        bool is_trivial();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        bool contains(NetworKit::node v);

        void check_consistency();
        void print(int depth = 0);
        void short_print();

        void delete_all_children();
        ~Blossom();
    };
    
    struct BacktrackInfo { 
        Blossom* blossom;  
        EdgeInfo edge; 
    };

    std::vector<bool> is_in_matching;
    std::vector<std::tuple<NetworKit::node, NetworKit::node, NetworKit::edgeweight>> graph_edges;
    std::set<Blossom*> blossoms;
    std::vector<Blossom*> trivial_blossom;
    std::vector<NetworKit::node> matched_vertex;

    void run_stage();
    virtual void initialize_stage() = 0;
    virtual void finish_stage() = 0;

    bool run_substage();
    virtual void initialize_substage() = 0;

    virtual bool has_useful_edges() = 0;
    virtual EdgeInfo get_useful_edge() = 0;

    bool consider_edge(EdgeInfo edge);
    virtual void label_odd(Blossom* b) = 0;
    virtual void label_even(Blossom* b) = 0;
    bool backtrack(Blossom* u, Blossom* v, EdgeInfo edge);
    bool backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path);
    
    void create_new_blossom(
            Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
    void cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut_2);
    virtual void handle_new_blossom(Blossom* b) = 0;

    void augment_path(
            Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
    void swap_edge_in_matching(NetworKit::edgeid edge);
    void swap_edges_on_even_path(
            Blossom* blossom, NetworKit::node u, NetworKit::node v);
    void augment_path_on_blossoms(
            Blossom* blossom, 
            std::vector<Blossom*>::reverse_iterator iter,
            NetworKit::node out_vertex);
    std::vector<Blossom*> blossoms_containing(NetworKit::node vertex, Blossom* until);
    static std::pair<std::list<std::pair<Blossom*, EdgeInfo>>,std::list<std::pair<Blossom*, EdgeInfo>>>
    split_subblossoms(std::list<std::pair<Blossom*, EdgeInfo>> sub_blossoms, Blossom* blossom);

    void adjust_dual_variables();
    
    virtual NetworKit::edgeweight calc_delta1() = 0;
    virtual NetworKit::edgeweight calc_delta2() = 0;
    virtual NetworKit::edgeweight calc_delta3() = 0;
    virtual NetworKit::edgeweight calc_delta4() = 0;
    virtual void adjust_by_delta(NetworKit::edgeweight delta) = 0;
    virtual void find_delta2_useful_edges() = 0;
    virtual void find_delta3_useful_edges() = 0;

    virtual void expand_odd_blossom(Blossom* blossom) = 0;
    virtual void expand_even_blossom(Blossom* blossom) = 0;

    bool is_exposed(Blossom* b);
    
    virtual Blossom* get_blossom(NetworKit::node vertex) = 0;

    static EdgeInfo reverse(const EdgeInfo& edge);
    static void print_backtrack(
            Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);

    virtual void check_consistency() = 0;
};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class EdmondsMaximumMatching final : public BlossomMaximumMatching {

public:
    EdmondsMaximumMatching(NetworKit::Graph &graph);

private:

    std::vector<Blossom*> current_blossom;
    std::vector<NetworKit::edgeweight> U;
    std::queue<EdgeInfo> useful_edges;

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void label_odd(Blossom* b) override;
    void label_even(Blossom* b) override;

    void handle_new_blossom(Blossom* b) override;
    
    NetworKit::edgeweight calc_delta1() override;
    NetworKit::edgeweight calc_delta2() override;
    NetworKit::edgeweight calc_delta3() override;
    NetworKit::edgeweight calc_delta4() override;
    void adjust_by_delta(NetworKit::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;

    void expand_odd_blossom(Blossom* blossom) override;
    void expand_even_blossom(Blossom* blossom) override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;

    NetworKit::edgeweight edge_dual_variable(NetworKit::edgeid edge);
    bool is_useful(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge);

};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class GabowMaximumMatching final : public BlossomMaximumMatching {

public:
    GabowMaximumMatching(NetworKit::Graph &graph);

private:
    class GabowBlossomData : public Blossom::BlossomData {
    public:
        GabowBlossomData(): best_edges(), best_edge(no_edge) {}
        std::unordered_map<Blossom*, EdgeInfo> best_edges;
        EdgeInfo best_edge;
    };
    GabowBlossomData* get_data(Blossom* b);
    
    std::queue<EdgeInfo> edge_queue;
    std::vector<NetworKit::edgeweight> U;
    std::vector<Blossom*> current_blossom;
    std::vector<EdgeInfo> best_edge;

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void label_odd(Blossom* b) override;
    void label_even(Blossom* b) override;

    void handle_new_blossom(Blossom* b) override;
    
    NetworKit::edgeweight calc_delta1() override;
    NetworKit::edgeweight calc_delta2() override;
    NetworKit::edgeweight calc_delta3() override;
    NetworKit::edgeweight calc_delta4() override;
    void adjust_by_delta(NetworKit::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;

    void expand_odd_blossom(Blossom* blossom) override;
    void expand_even_blossom(Blossom* blossom) override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;

    void calc_best_edges(Blossom* b);

    NetworKit::edgeweight edge_slack(NetworKit::edgeid edge);
};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class MicaliGabowMaximumMatching final : public BlossomMaximumMatching {

public:
    MicaliGabowMaximumMatching(NetworKit::Graph &graph);

private:
    using BlossomNodeList = ConcatenableQueue<Blossom*, NetworKit::node, NetworKit::node>;
    using EvenEdgeGroup = PriorityQueue2<NetworKit::edgeid, NetworKit::edgeweight>::Group*;
    class MicaliGabowBlossomData : public Blossom::BlossomData {
    public:
        // Contains all nodes in blossom order
        BlossomNodeList nodes;
        // Group corresponding to blossom in queue even_edges
        EvenEdgeGroup even_edges;

        MicaliGabowBlossomData(BlossomNodeList&& nodes, EvenEdgeGroup even_edges): 
            nodes(std::move(nodes)), even_edges(even_edges) {}
    };
    MicaliGabowBlossomData* get_data(Blossom* b);
    // References to nodes in concatenable queues of blossoms
    std::vector<BlossomNodeList::ElementRef> nodes_refs;

    // Queue of useful edges
    std::queue<EdgeInfo> edge_queue;

    // Used to maintain dual variables for vertices
    PriorityQueue1<NetworKit::node, NetworKit::edgeweight> Ueven;
    PriorityQueue1<NetworKit::node, NetworKit::edgeweight> Uodd;
    std::vector<NetworKit::edgeweight> Ufree;
    NetworKit::edgeweight U(NetworKit::node v);

    // Used to maintain dual variables for blossoms
    PriorityQueue1<NetworKit::node, NetworKit::edgeweight> Zeven;
    PriorityQueue1<NetworKit::node, NetworKit::edgeweight> Zodd;
    NetworKit::edgeweight blossom_dual(Blossom* b);

    // Used to maintain slack of edges between S-vertices
    // Needed to calculate delta_3
    // Maintains pi_ij / 2
    PriorityQueue1<NetworKit::edgeid, NetworKit::edgeweight> good_edges;
    void clear_not_good_edges();
    bool is_good(NetworKit::edgeid edge);
    NetworKit::edgeweight edge_slack(NetworKit::edgeid edge);

    // Used to maintain edges from S-vertices to T/free-vertices
    // Neede to calculate delta_2
    // Maintains pi_ij for 
    PriorityQueue2<NetworKit::edgeid, NetworKit::edgeweight> even_edges;
    NetworKit::edgeid dummy_edge_id(NetworKit::node node);

    // Scan all edges from newly even blossom and update good_edges and even_edges
    void scan_edges(Blossom* b);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void label_odd(Blossom* b) override;
    void label_even(Blossom* b) override;

    void handle_new_blossom(Blossom* b) override;
    
    NetworKit::edgeweight calc_delta1() override;
    NetworKit::edgeweight calc_delta2() override;
    NetworKit::edgeweight calc_delta3() override;
    NetworKit::edgeweight calc_delta4() override;
    void adjust_by_delta(NetworKit::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;

    void expand_odd_blossom(Blossom* blossom) override;
    void expand_even_blossom(Blossom* blossom) override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;
};

} /* namespace Koala */
