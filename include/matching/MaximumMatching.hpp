#pragma once

#include <map>
#include <set>
#include <list>
#include <deque>
#include <string>
#include <optional>
#include <functional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <matching/PriorityQueues.hpp>

#define DEBUG_LOGGING 0

namespace Koala {

/**
 * @ingroup matching
 * The base class for the maximum matching algorithms.
 *
 */
class MaximumMatching : public NetworKit::Algorithm {

public:
    using edgeweight = NetworKit::edgeweight;
    using intedgeweight = int;

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

        bool operator==(const EdgeInfo& other) const { return other.id == id; }
    };

    static constexpr EdgeInfo no_edge { NetworKit::none, NetworKit::none, NetworKit::none };
    
    static EdgeInfo reverse(const EdgeInfo& edge);

    std::string edge_to_string(const EdgeInfo& e) { 
        return (e == no_edge) ? 
            "none" : 
            "(" + std::to_string(e.u) + ", " + std::to_string(e.v) + ")"; 
    }

    enum BlossomLabel { odd, even, free };

    class Blossom {
    public:
        class BlossomData {};
        Blossom* parent;
        
        // Base of the blossom at the moment of creation
        NetworKit::node initial_base;
        
        // Current base of the blossom
        NetworKit::node base;
        
        // Last node in the blossom order
        NetworKit::node last_node;
        
        std::list<Blossom*> base_blossoms;
        
        // Subblossoms of the blossom
        std::list<std::pair<Blossom*, EdgeInfo>> subblossoms;
        
        // Label and backtrack edge set during a search
        BlossomLabel label;
        EdgeInfo backtrack_edge;
        
        // Visited flag used during backtracking
        bool visited;
        
        // Dual weight corresponding to the blossom
        // To be used by the specific algorithms, might not be correct depending on implementation
        MaximumMatching::edgeweight z;

        // Iterator pointing to the blossom in the 
        std::list<Blossom*>::iterator list_it;

        // Pointer to algorithm-specific data corresponding to a blossom
        BlossomData* data;

        bool is_trivial();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        bool contains(NetworKit::node v);

        void check_consistency();
        void print(int depth = 0);
        void short_print();
        void nodes_print();

        void delete_all_children();
        ~Blossom();
    };
    
    struct BacktrackInfo { 
        Blossom* blossom;  
        EdgeInfo edge; 
    };

    // Store all graph edges so they can be accessed using the id
    std::vector<std::tuple<NetworKit::node, NetworKit::node, MaximumMatching::edgeweight>> graph_edges;
    
    // Store a list of all proper blossoms
    std::list<Blossom*> blossoms;
    
    // For each edge store wether it's in the current matching
    std::vector<bool> is_in_matching;
    
    // For each vertex store the vertex it's matched to and the id of the edge by which they're matched
    std::vector<NetworKit::node> matched_vertex;
    std::vector<NetworKit::edgeid> matched_edge;

    // For each vertex store it's corresponding trivial blossom
    std::vector<Blossom*> trivial_blossom;

    void run_stage();
    virtual void initialize_stage() = 0;
    virtual void finish_stage() = 0;

    bool run_substage();
    virtual void initialize_substage() = 0;

    void expand_final_blossom(Blossom* b);

    virtual bool has_useful_edges() = 0;
    virtual EdgeInfo get_useful_edge() = 0;

    bool consider_edge(EdgeInfo edge);

    virtual void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) = 0;

    bool backtrack(Blossom* u, Blossom* v, EdgeInfo edge);
    bool backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path);
    
    void add_blossom(Blossom* b);
    void remove_blossom(Blossom* b);
    void create_new_blossom(Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
    void cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut_2);
    virtual void handle_new_blossom(Blossom* b) = 0;

    void augment_path(Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);

    void swap_edges_on_even_path(Blossom* blossom, NetworKit::node u, NetworKit::node v);
    void swap_edges_on_even_path(Blossom* blossom, NetworKit::node out_vertex, std::list<Blossom*>&& base_blossoms);

    void lazy_augment_path_in_blossom(Blossom* blossom);

    std::pair<std::list<std::pair<Blossom*, EdgeInfo>>,std::list<std::pair<Blossom*, EdgeInfo>>>
    split_subblossoms(std::list<std::pair<Blossom*, EdgeInfo>> subblossoms, Blossom* blossom);
    void swap_edge_in_matching(NetworKit::edgeid edge);
    void check_edge_in_matching(NetworKit::edgeid edge);
    std::list<Blossom*> blossoms_containing(NetworKit::node u, Blossom* until);

    virtual void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) = 0;

    void adjust_dual_variables();
    
    virtual MaximumMatching::edgeweight calc_delta1() = 0;
    virtual MaximumMatching::edgeweight calc_delta2() = 0;
    virtual MaximumMatching::edgeweight calc_delta3() = 0;
    virtual MaximumMatching::edgeweight calc_delta4() = 0;
    
    virtual void adjust_by_delta(MaximumMatching::edgeweight delta) = 0;
    virtual void find_delta2_useful_edges() = 0;
    virtual void find_delta3_useful_edges() = 0;
    virtual std::vector<Blossom*> get_odd_blossoms_to_expand() = 0;

    void expand_odd_blossom(Blossom* blossom);
    void expand_even_blossom(Blossom* blossom);
    
    virtual void handle_odd_blossom_expansion(Blossom* blossom) = 0;
    virtual void handle_even_blossom_expansion(Blossom* blossom) = 0;

    bool is_exposed(Blossom* b);
    
    virtual Blossom* get_blossom(NetworKit::node vertex) = 0;
    
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

    // For each vertex store it's current blossom and the value of corresponding dual variable
    std::vector<Blossom*> current_blossom;
    std::vector<MaximumMatching::edgeweight> U;
    
    // Queue of tight edges
    std::queue<EdgeInfo> edge_queue;

    void scan_edges(Blossom* b);

    MaximumMatching::edgeweight slack(NetworKit::edgeid edge);
    bool is_tight(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;

    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;
    
    MaximumMatching::edgeweight calc_delta1() override;
    MaximumMatching::edgeweight calc_delta2() override;
    MaximumMatching::edgeweight calc_delta3() override;
    MaximumMatching::edgeweight calc_delta4() override;
    void adjust_by_delta(MaximumMatching::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;
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

        // For an even blossom store a list of edges (x_i, y_i) such that
        //    1. y_1 < y_2 < ... < y_n
        //    2. x_i belongs to the blossom
        //    3. y_i belongs to a different even blossom
        //    4. x_i became even after y_i
        std::list<EdgeInfo> best_edges;

        // Edge in the best_edges list with smallest slack
        EdgeInfo best_edge;
    };
    GabowBlossomData* get_data(Blossom* b);

    // For each vertex store a sorted list of it's neighbours
    // This is needed to efficiently create lists of best edges
    std::vector<std::vector<std::pair<NetworKit::node, NetworKit::edgeid>>> sorted_neighbours;
    
    // Queue of tight edges
    std::queue<EdgeInfo> edge_queue;

    // For each vertex store it's current blossom and the value of it's corresponding dual variable
    std::vector<Blossom*> current_blossom;
    std::vector<MaximumMatching::edgeweight> U;

    // For each non even vertex store an edge connecting it to an even vertex with minimum slack
    std::vector<EdgeInfo> best_edge;

    void scan_edges(Blossom* b);

    MaximumMatching::edgeweight edge_slack(NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;
    
    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;

    MaximumMatching::edgeweight calc_delta1() override;
    MaximumMatching::edgeweight calc_delta2() override;
    MaximumMatching::edgeweight calc_delta3() override;
    MaximumMatching::edgeweight calc_delta4() override;
    void adjust_by_delta(MaximumMatching::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;
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
    using EvenEdgeGroup = PriorityQueue2<NetworKit::edgeid, MaximumMatching::edgeweight>::Group*;
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
    PriorityQueue1<NetworKit::node, MaximumMatching::edgeweight> Ueven;
    PriorityQueue1<NetworKit::node, MaximumMatching::edgeweight> Uodd;
    std::vector<MaximumMatching::edgeweight> Ufree;
    MaximumMatching::edgeweight U(NetworKit::node v);

    // Used to maintain dual variables for blossoms
    PriorityQueue1<NetworKit::node, MaximumMatching::edgeweight> Zeven;
    PriorityQueue1<NetworKit::node, MaximumMatching::edgeweight> Zodd;
    MaximumMatching::edgeweight blossom_dual(Blossom* b);

    // Used to maintain slack of edges between S-vertices
    // Needed to calculate delta_3
    // Maintains pi_ij / 2
    PriorityQueue1<NetworKit::edgeid, MaximumMatching::edgeweight> good_edges;
    void clear_not_good_edges();
    bool is_good(NetworKit::edgeid edge);

    // Used to maintain edges from S-vertices to T/free-vertices
    // Needed to calculate delta_2
    // Maintains pi_ij for 
    PriorityQueue2<NetworKit::edgeid, MaximumMatching::edgeweight> even_edges;
    NetworKit::edgeid dummy_edge_id(NetworKit::node node);

    // Scan all edges from newly even blossom and update good_edges and even_edges
    void scan_edges(Blossom* b);

    MaximumMatching::edgeweight edge_slack(NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    EdgeInfo get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;

    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;
    
    MaximumMatching::edgeweight calc_delta1() override;
    MaximumMatching::edgeweight calc_delta2() override;
    MaximumMatching::edgeweight calc_delta3() override;
    MaximumMatching::edgeweight calc_delta4() override;
    void adjust_by_delta(MaximumMatching::edgeweight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;
    
    Blossom* get_blossom(NetworKit::node vertex) override;

    void check_consistency() override;
};

/**
 * @ingroup matching
 * The class for Gabow's scaling maximum matching algorithm.
 */
class GabowScalingMatching : public MaximumMatching {
public:
    GabowScalingMatching(NetworKit::Graph &graph);

    /**
     * Execute the Gabow scaling maximum matching algorithm.
     */
    void run();

private:
    struct EdgeInfo { 
        NetworKit::node u, v; 
        NetworKit::edgeid id; 

        bool operator==(const EdgeInfo& other) const;
        bool operator!=(const EdgeInfo& other) const;
        bool operator<(const EdgeInfo& other) const;
    };

    friend std::ostream& operator<<(std::ostream &out, const EdgeInfo& edge);
    static constexpr EdgeInfo no_edge {NetworKit::none, NetworKit::none, NetworKit::none};

    enum BlossomLabel { even, odd, free };

    struct Blossom {
        NetworKit::node base;
        MaximumMatching::intedgeweight z0, t_root, t_odd, t_even, Delta;
        Blossom* parent;
        std::list<std::pair<Blossom*, EdgeInfo>> subblossoms;
        BlossomLabel label;
        EdgeInfo backtrack_edge;
        bool visited;
        std::list<Blossom*>::iterator shell_blossoms_it;
        SplitFindMin<NetworKit::node, Blossom*, MaximumMatching::intedgeweight, EdgeInfo>::List* list;

        bool is_trivial();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        std::list<NetworKit::node> node_list();
        
        void short_print();
        void nodes_print();
    };

    struct Shell {
        // The number of all nodes in the shell
        int size;

        // Dual weight of the shell
        MaximumMatching::intedgeweight z;

        // TODO - check if these can be moved outside
        
        // Time the shell has become the active shell in the search
        MaximumMatching::intedgeweight t_active;

        // Time the shell has become the inner shell in the search
        MaximumMatching::intedgeweight t_inner;

        std::list<Shell*> children;

        // List of nodes solely in this shell
        std::list<NetworKit::node> nodes;

        // Parent and child in heavy path decomposition
        Shell* heavy_path_parent;
        Shell* heavy_child;

        // Blossoms contained in the shell
        std::list<Blossom*> shell_blossoms;

        // Status of the shell - wether it has been searched or dissolved
        bool dissolved, searched;

        // Index of the shell in the path
        int shell_index;

        bool is_heavy_path_root();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        void for_blossoms(const std::function<void(Shell*)>& handle);
        void short_print();
    };

    NetworKit::Graph reducedGraph;
    std::vector<std::pair<NetworKit::node, NetworKit::node>> graph_edges; 

    std::tuple<std::vector<NetworKit::node>, std::vector<MaximumMatching::intedgeweight>, Shell*>
    scale(const std::vector<MaximumMatching::intedgeweight>& w);

    void match(Shell* T);

    std::vector<MaximumMatching::intedgeweight> current_w;
    std::vector<MaximumMatching::intedgeweight> current_y;
    std::vector<Blossom*> trivial_blossom;
    std::vector<NetworKit::node> matched_vertex;
    std::vector<NetworKit::edgeid> matched_edge;
    std::vector<bool> edge_in_matching;
    int edges_in_matching;

    void delete_blossom(Blossom* B, Shell* S);
    void add_blossom(Blossom* B, Shell* S);
    void expand_blossom(Blossom* B, Shell* S);
    void swap_edge_in_matching(NetworKit::edgeid edge);
    void set_edge_in_matching(NetworKit::edgeid edge);
    void remove_edge_from_matching(NetworKit::edgeid edge);
    void check_edge_in_matching(NetworKit::edgeid edge);

    void create_trivial_blossoms(Shell* T);
    Shell* turn_current_blossoms_into_old(const std::list<Blossom*>& subblossoms);
    void heavy_path_decomposition(Shell* T, MaximumMatching::intedgeweight outer_dual);
    void change_blossom_base(Blossom* B, NetworKit::node new_base, EdgeInfo edge);
    void swap_edges_on_even_path(Blossom* B, NetworKit::node u, NetworKit::node v);
    void swap_edges_on_even_path(Blossom* B, NetworKit::node out_vertex, std::list<Blossom*>&& out_blossoms);
    
    void path(Shell* B, MaximumMatching::intedgeweight outer_dual);

    int free_nodes_in_shell(Shell* B);
    void enumerate_shells();
    void augmentPaths();
    int current_slack(EdgeInfo e);
    
    std::vector<NetworKit::node> actual_to_contracted;

    struct Event {
        enum Type { grow, blossom, dissolve, dissolveShell };
        Type type;
        union {
            struct { NetworKit::node v; EdgeInfo e; };
            EdgeInfo uv;
            Blossom* B;
            Shell* S;
        } args;

        static Event make_grow(NetworKit::node v, EdgeInfo e) { 
            Event event;
            event.type = grow;
            event.args.v = v; event.args.e = e; 
            return event;
        }

        static Event make_blossom(EdgeInfo uv) { 
            Event event;
            event.type = blossom;
            event.args.uv = uv; 
            return event;
        }

        static Event make_dissolve(Blossom* B) { 
            Event event;
            event.type = dissolve;
            event.args.B = B;
            return event;
        }

        static Event make_dissolveShell(Shell* S) { 
            Event event;
            event.type = dissolveShell;
            event.args.S = S; 
            return event;
        }
    };

    friend std::ostream& operator<<(std::ostream &out, const Event& event);

    void shell_search(Shell* B);
    void init_shell_search();
    
    void finish_shell_search();
    void delete_lists(Blossom* B);

    bool search_done;
    MaximumMatching::intedgeweight t_undissolvable;
    MaximumMatching::intedgeweight shell_Delta;
    MaximumMatching::intedgeweight outer_blossom_dual;
    MaximumMatching::intedgeweight active_shell_initial_dual;
    Shell* old_root;
    Shell* path_root;
    Shell* highest_undissolved;
    Shell* active_shell;
    Shell* starting_shell;
    std::vector<std::pair<MaximumMatching::intedgeweight, Shell*>> shells;
    std::vector<Shell*> vertex_path;
    std::vector<MaximumMatching::intedgeweight> current_shell_duals;
    std::vector<MaximumMatching::intedgeweight> y0;
    std::vector<MaximumMatching::intedgeweight> t_shell;
    std::vector<MaximumMatching::intedgeweight> Delta;
    std::vector<Blossom*> current_blossom;
    std::vector<Shell*> current_shell;
    std::vector<Shell*> search_shell;
    FenwickTree shell_distribution;
    ArrayPriorityQueue<Event> event_queue;
    UnionFind<NetworKit::node, Blossom*> union_find;
    SplitFindMin<NetworKit::node, Blossom*, MaximumMatching::intedgeweight, EdgeInfo> split_find_min;

    void grow(NetworKit::node v, EdgeInfo e);
    void schedule(NetworKit::node v);
    void dissolve(Blossom* B);
    void blossom(EdgeInfo edge);
    void dissolveShell(Shell* S);

    struct BacktrackInfo { 
        Blossom* blossom;  
        EdgeInfo edge; 
    };

    bool backtrack(Blossom* u, Blossom* v, EdgeInfo edge);
    bool backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path);
    
    void create_new_blossom(Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
    
    std::list<Blossom*> blossoms_containing(NetworKit::node u, Blossom* until);
    void cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut_2);

    std::pair<std::list<std::pair<Blossom*, EdgeInfo>>, std::list<std::pair<Blossom*, EdgeInfo>>>
    split_subblossoms(std::list<std::pair<Blossom*, EdgeInfo>> subblossoms, Blossom* blossom);

    MaximumMatching::intedgeweight y(NetworKit::node v);
    MaximumMatching::intedgeweight z(Blossom* B);
    MaximumMatching::intedgeweight shell_z();
    MaximumMatching::intedgeweight slack(EdgeInfo e);
    bool is_exposed(Blossom* b);
    bool is_even(NetworKit::node v);
    bool is_odd(NetworKit::node v);
    Blossom* get_blossom(NetworKit::node v);
    void add_distribution(Shell* S, MaximumMatching::intedgeweight distribution);
    MaximumMatching::intedgeweight distribution_so_far(int shell_index);
    bool matching_is_perfect();

    void print_current_state();
    void print_search_state();
    void print_vertex_state(NetworKit::node v);
    void check_consistency();
    void check_consistency(Blossom* B);

    static std::string label_to_str(BlossomLabel label);
    static EdgeInfo reverse(const EdgeInfo& edge);
    static void print_backtrack(
            Blossom* u, Blossom* v, EdgeInfo edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
};

/**
 * @ingroup matching
 * The base class for the maximum cardinality matching algorithms.
 *
 */
class MaximumCardinalityMatching : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the maximum matching algorithm.
     *
     * @param graph The input graph.
     */
    MaximumCardinalityMatching(NetworKit::Graph &graph);

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
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class MicaliVaziraniMatching final : public MaximumCardinalityMatching {

public:
    MicaliVaziraniMatching(NetworKit::Graph &graph);
    MicaliVaziraniMatching(
        NetworKit::Graph &graph, const std::vector<NetworKit::node>& initial_matching);

    /**
     * Execute the Micali Vazirani maximum cardinality matching algorithm.
     */
    void run();

private:
    struct Bloom {
        NetworKit::node base;
        int green_color;
        int red_color;
        NetworKit::node green_peak;
        NetworKit::node green_root;
        NetworKit::node red_peak;
        NetworKit::node red_root;
    };

    static constexpr int inf_level = 1e9;

    struct VertexData {
        NetworKit::node match;
        NetworKit::node parent;
        int even_level;
        int odd_level;
        Bloom* bloom;
        int color;
        std::vector<NetworKit::node> predecessors;
        size_t pred_it;
        int pred_count;
        std::vector<NetworKit::node> successors;
        std::vector<std::pair<NetworKit::node, NetworKit::node>> children;
        bool erased;
    };

    static constexpr int no_color = 0;
    int color_counter;

    struct EdgeData {
        enum Type { none, prop, bridge };
        NetworKit::node u, v;
        Type type;
    };

    std::vector<VertexData> V;
    std::vector<EdgeData> E;
    std::vector<std::vector<NetworKit::node>> candidates;
    std::vector<std::vector<NetworKit::edgeid>> bridges;

    bool augmentation_happened;
    bool bloom_found;
    int iter, max_iter;

    std::vector<Bloom*> current_blooms;
    UnionFind<NetworKit::node, NetworKit::node> bloom_bases;
    std::vector<NetworKit::node> bridge_support;
    std::vector<NetworKit::node> erase_queue;

    void reset();
    void search();
    void clear_blooms();

    void bloss_aug(NetworKit::node s, NetworKit::node t);
    void red_dfs_step(NetworKit::node& v_R, int red_color, NetworKit::node& v_G, 
                  NetworKit::node r, NetworKit::node& barrier);
    void green_dfs_step(NetworKit::node& v_G, int green_color, NetworKit::node& v_R, 
                  NetworKit::node& barrier);

    void erase(std::vector<NetworKit::node>& Y);
    
    void find_path(NetworKit::node x, NetworKit::node y);
    bool open(NetworKit::node cur, NetworKit::node bcur, NetworKit::node b);
    
    NetworKit::node base_star(Bloom* bloom);
    NetworKit::node base_star(NetworKit::node vertex);
    NetworKit::node base(NetworKit::node vertex);

    void flip_edge(NetworKit::node u, NetworKit::node v);
    void set_level(NetworKit::node vertex, int level);
    bool exposed(NetworKit::node vertex);
    int min_level(NetworKit::node vertex);
    int max_level(NetworKit::node vertex);
    int tenacity(NetworKit::node u, NetworKit::node v);
    bool outer(NetworKit::node vertex);
    bool inner(NetworKit::node vertex);

    std::tuple<NetworKit::node, NetworKit::node, NetworKit::node, NetworKit::node> 
    get_bridge(NetworKit::node vertex);
};

} /* namespace Koala */
