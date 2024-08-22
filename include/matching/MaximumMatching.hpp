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

#include <matching/structures/ArrayPriorityQueue.hpp>
#include <matching/structures/ConcatenableQueue.hpp>
#include <matching/structures/FenwickTree.hpp>
#include <matching/structures/PriorityQueues.hpp>
#include <matching/structures/SplitFindMin.hpp>
#include <matching/structures/UnionFind.hpp>

namespace Koala {

/**
 * @ingroup matching
 * The base class for the maximum weight matching algorithms.
 *
 */
class MaximumWeightMatching :  public NetworKit::Algorithm {
 public:
    // This is double, which works for non-scaling algorithms
    using weight = NetworKit::edgeweight;
    using intweight = int;

    /**
     * Given an input graph, set up the maximum matching algorithm.
     *
     * @param graph The input graph.
     */
    explicit MaximumWeightMatching(NetworKit::Graph &graph, bool perfect = false) :
        graph(graph), perfect(perfect) {}

    /**
     * Return the matching found by the algorithm.
     *
     * @return a matching between nodes.
     */
    const std::map<NetworKit::node, NetworKit::node>& getMatching() const {
        assureFinished();
        return matching;
    }

 protected:
    static constexpr weight infinite_weight = std::numeric_limits<weight>::max();

    NetworKit::Graph graph;
    bool perfect;

    std::map<NetworKit::node, NetworKit::node> matching;
};

/**
 * Abstract class for implementation's of Edmond's Blossom algorithm.
 * Implements essential steps in the algorithm while deferring some calculations to it's child
 * class which implement different variants of the algorithm.
*/
class BlossomMaximumMatching :  public MaximumWeightMatching {
 public:
    explicit BlossomMaximumMatching(NetworKit::Graph &graph, bool perfect);

    /**
     * Execute the Edmonds maximum matching algorithm.
     */
    void run();

 protected:
    virtual ~BlossomMaximumMatching();

    struct Edge {
        NetworKit::node u, v;
        NetworKit::edgeid id;

        bool operator==(const Edge& other) const { return other.id == id; }
    };

    static constexpr Edge no_edge { NetworKit::none, NetworKit::none, NetworKit::none };

    static Edge reverse(const Edge& edge);
    static std::string edge_to_string(const Edge& e);

    enum Label { odd, even, free };

    class Blossom {
     public:
        using SubblossomList = std::list<std::pair<Blossom*, Edge>>;

        class BlossomData { public: virtual ~BlossomData() {} };

        static Blossom* trivial(NetworKit::node vertex);
        static Blossom* nontrivial(const SubblossomList& subblossoms);

        Blossom* parent;

        // Base of the blossom at the moment of creation
        NetworKit::node initial_base;

        // Current base of the blossom
        NetworKit::node base;

        // Last node in the blossom order
        NetworKit::node last_node;

        // List of the descendants containg the base, used for lazy augmentation
        std::list<Blossom*> base_blossoms;

        // Subblossoms of the blossom
        SubblossomList subblossoms;

        // List of nodes in a root blossom
        std::list<NetworKit::node> nodes;

        // Label and backtrack edge set during a search
        Label label;
        Edge backtrack_edge;

        // Visited flag used during backtracking
        bool visited;

        // Dual weight corresponding to the blossom
        // To be used by the specific algorithms, might not be correct depending on implementation
        MaximumWeightMatching::weight z;

        // Iterator pointing to the blossom's position in the root blossom list
        std::list<Blossom*>::iterator list_it;

        // Pointer to algorithm-specific data corresponding to a blossom
        BlossomData* data;

        bool is_trivial();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        bool contains(NetworKit::node v);

        void delete_all_children();
        ~Blossom();
    };

    struct BacktrackInfo {
        Blossom* blossom;
        Edge edge;
    };

    bool finished;

    // Maximum edge weight in the graph
    MaximumWeightMatching::weight max_weight;

    // Store all graph edges so they can be accessed using the id in constant time
    std::vector<std::tuple<NetworKit::node, NetworKit::node, MaximumWeightMatching::weight>>
    graph_edges;

    // Store a list of all proper blossoms
    std::list<Blossom*> blossoms;

    // Store iterators to node positions in blossom node lists
    std::vector<std::list<NetworKit::node>::iterator> node_iter;

    // For each edge store wether it's in the current matching
    std::vector<bool> is_in_matching;
    NetworKit::count edges_in_matching;

    // For each vertex store the vertex it's match and the id of the edge by which they're matched
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
    virtual Edge get_useful_edge() = 0;

    bool consider_edge(Edge edge);

    virtual void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) = 0;

    bool backtrack(Blossom* u, Blossom* v, Edge edge);
    bool backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path);

    void add_blossom(Blossom* b);
    void remove_blossom(Blossom* b);
    void create_new_blossom(Blossom* u, Blossom* v, Edge edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
    void cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut_2);
    virtual void handle_new_blossom(Blossom* b) = 0;

    void augment_path(Blossom* u, Blossom* v, Edge edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);

    void swap_edges_on_even_path(Blossom* blossom, NetworKit::node u, NetworKit::node v);
    void swap_edges_on_even_path(Blossom* blossom, NetworKit::node out_vertex,
        std::list<Blossom*>&& base_blossoms);

    void lazy_augment_path_in_blossom(Blossom* blossom);

    std::pair<Blossom::SubblossomList, Blossom::SubblossomList>
    split_subblossoms(Blossom::SubblossomList subblossoms, Blossom* blossom);
    void swap_edge_in_matching(NetworKit::edgeid edge);
    void check_edge_in_matching(NetworKit::edgeid edge);
    std::list<Blossom*> blossoms_containing(NetworKit::node u, Blossom* until);

    virtual void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) = 0;

    void adjust_dual_variables();

    virtual MaximumWeightMatching::weight calc_delta1() = 0;
    virtual MaximumWeightMatching::weight calc_delta2() = 0;
    virtual MaximumWeightMatching::weight calc_delta3() = 0;
    virtual MaximumWeightMatching::weight calc_delta4() = 0;

    virtual void adjust_by_delta(MaximumWeightMatching::weight delta) = 0;
    virtual void find_delta2_useful_edges() = 0;
    virtual void find_delta3_useful_edges() = 0;
    virtual std::vector<Blossom*> get_odd_blossoms_to_expand() = 0;

    void expand_odd_blossom(Blossom* blossom);
    void expand_even_blossom(Blossom* blossom);

    virtual void handle_odd_blossom_expansion(Blossom* blossom) = 0;
    virtual void handle_even_blossom_expansion(Blossom* blossom) = 0;

    bool is_exposed(Blossom* b);

    bool is_matching_perfect();

    virtual Blossom* get_blossom(NetworKit::node vertex) = 0;
};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class EdmondsMaximumMatching final :  public BlossomMaximumMatching {
 public:
    explicit EdmondsMaximumMatching(NetworKit::Graph &graph, bool perfect = true);

 private:
    // For each vertex store it's current blossom and the value of corresponding dual variable
    std::vector<Blossom*> current_blossom;
    std::vector<MaximumWeightMatching::weight> y;

    // Queue of tight edges
    std::queue<Edge> edge_queue;

    void scan_edges(Blossom* b);

    MaximumWeightMatching::weight slack(NetworKit::edgeid edge);
    bool is_tight(NetworKit::node u, NetworKit::node v, NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    Edge get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;

    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;

    MaximumWeightMatching::weight calc_delta1() override;
    MaximumWeightMatching::weight calc_delta2() override;
    MaximumWeightMatching::weight calc_delta3() override;
    MaximumWeightMatching::weight calc_delta4() override;
    void adjust_by_delta(MaximumWeightMatching::weight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;

    Blossom* get_blossom(NetworKit::node vertex) override;
};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class GabowMaximumMatching final :  public BlossomMaximumMatching {
 public:
    explicit GabowMaximumMatching(NetworKit::Graph &graph, bool perfect = true);

 private:
    class GabowBlossomData :  public Blossom::BlossomData {
     public:
        GabowBlossomData(): best_edges(), best_edge(no_edge) {}

        // For an even blossom store a list of edges (x_i, y_i) such that
        //    1. y_1 < y_2 < ... < y_n
        //    2. x_i belongs to the blossom
        //    3. y_i belongs to a different even blossom
        //    4. x_i became even after y_i
        std::list<Edge> best_edges;

        // Edge in the best_edges list with smallest slack
        Edge best_edge;

        ~GabowBlossomData() override = default;
    };
    GabowBlossomData* get_data(Blossom* b);

    // Queue of tight edges
    std::queue<Edge> edge_queue;

    // For each vertex store it's current blossom and the value of it's corresponding dual variable
    std::vector<Blossom*> current_blossom;
    std::vector<MaximumWeightMatching::weight> y;

    // For each non even vertex store an edge connecting it to an even vertex with minimum slack
    std::vector<Edge> best_edge;

    void scan_edges(Blossom* b);

    MaximumWeightMatching::weight slack(NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    Edge get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;

    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;

    MaximumWeightMatching::weight calc_delta1() override;
    MaximumWeightMatching::weight calc_delta2() override;
    MaximumWeightMatching::weight calc_delta3() override;
    MaximumWeightMatching::weight calc_delta4() override;
    void adjust_by_delta(MaximumWeightMatching::weight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;

    Blossom* get_blossom(NetworKit::node vertex) override;
};

/**
 * @ingroup matching
 * The class for the Edmonds maximum matching algorithm.
 */
class GalilMicaliGabowMaximumMatching final :  public BlossomMaximumMatching {
 public:
    explicit GalilMicaliGabowMaximumMatching(NetworKit::Graph &graph, bool perfect = true);

 private:
    using BlossomNodeList = ConcatenableQueue<NetworKit::node, NetworKit::node, Blossom*>;
    using EvenEdgeQueue = PriorityQueue2<NetworKit::node, MaximumWeightMatching::weight,
        NetworKit::edgeid>;
    using EvenEdgeGroup = EvenEdgeQueue::Group*;

    class GalilMicaliGabowBlossomData : public Blossom::BlossomData {
     public:
        // Contains all nodes in blossom order
        BlossomNodeList nodes;
        // Group corresponding to blossom in queue even_edges
        EvenEdgeGroup even_edge_group;

        GalilMicaliGabowBlossomData(BlossomNodeList&& nodes, EvenEdgeGroup even_edge_group):
            nodes(std::move(nodes)), even_edge_group(even_edge_group) {}

        ~GalilMicaliGabowBlossomData() override = default;
    };
    GalilMicaliGabowBlossomData* get_data(Blossom* b);
    // References to nodes in concatenable queues of blossoms
    std::vector<BlossomNodeList::handle_type> nodes_refs;

    // Queue of useful edges
    std::queue<Edge> edge_queue;

    // Used to maintain dual variables for vertices
    PriorityQueue1<NetworKit::node, MaximumWeightMatching::weight, NetworKit::node> y_even;
    PriorityQueue1<NetworKit::node, MaximumWeightMatching::weight, NetworKit::node> y_odd;
    std::vector<MaximumWeightMatching::weight> y_free;
    MaximumWeightMatching::weight y(NetworKit::node v);

    // Used to maintain dual variables for blossoms
    PriorityQueue1<NetworKit::node, MaximumWeightMatching::weight, Blossom*> z_even;
    PriorityQueue1<NetworKit::node, MaximumWeightMatching::weight, Blossom*> z_odd;
    MaximumWeightMatching::weight z(Blossom* b);

    // Used to maintain slack of edges between even vertices in different blossoms
    // Needed to calculate delta_3
    // Maintains slack / 2 for those edges
    // During the execution some edges might become contained in the same even blossom and have to
    // be removed. This is done lazily when checking for a good edge with smallest slack.
    PriorityQueue1<NetworKit::edgeid, MaximumWeightMatching::weight, NetworKit::edgeid> good_edges;
    void clear_not_good_edges();
    bool is_good(NetworKit::edgeid edge);

    // Used to maintain edges from even to non-even vertices
    // Needed to calculate delta_2
    // Maintains slack for those edges
    EvenEdgeQueue even_edges;

    // Scan all edges from newly even blossom and update good_edges and even_edges
    void scan_edges(Blossom* b);

    MaximumWeightMatching::weight slack(NetworKit::edgeid edge);

    void initialize_stage() override;
    void finish_stage() override;
    void initialize_substage() override;

    bool has_useful_edges() override;
    Edge get_useful_edge() override;

    void handle_grow(Blossom* odd_blossom, Blossom* even_blossom) override;

    void handle_new_blossom(Blossom* b) override;

    void handle_subblossom_shift(Blossom* blossom, Blossom* subblossom) override;
    void handle_odd_blossom_expansion(Blossom* blossom) override;
    void handle_even_blossom_expansion(Blossom* blossom) override;

    MaximumWeightMatching::weight calc_delta1() override;
    MaximumWeightMatching::weight calc_delta2() override;
    MaximumWeightMatching::weight calc_delta3() override;
    MaximumWeightMatching::weight calc_delta4() override;
    void adjust_by_delta(MaximumWeightMatching::weight delta) override;

    void find_delta2_useful_edges() override;
    void find_delta3_useful_edges() override;
    std::vector<Blossom*> get_odd_blossoms_to_expand() override;

    Blossom* get_blossom(NetworKit::node vertex) override;
};

/**
 * @ingroup matching
 * The class for Gabow's scaling maximum matching algorithm.
 */
class GabowScalingMatching :  public MaximumWeightMatching {
 public:
    explicit GabowScalingMatching(NetworKit::Graph &graph, bool perfect = false);

    /**
     * Execute the Gabow scaling maximum matching algorithm.
     */
    void run();

 private:
    struct Edge {
        NetworKit::node u, v;
        NetworKit::edgeid id;

        bool operator==(const Edge& other) const;
        bool operator!=(const Edge& other) const;
        bool operator<(const Edge& other) const;
    };

    friend std::ostream& operator<<(std::ostream &out, const Edge& edge);
    static constexpr Edge no_edge {NetworKit::none, NetworKit::none, NetworKit::none};

    enum Label { even, odd, free };

    struct Blossom {
        // The blossom's base
        NetworKit::node base;

        // Counters used in search:
        //  z0     - initial value of the dual variable
        //  t_root - time the blossom became a root blossom
        //  t_odd  - time the blossom became odd
        //  t_even - time the blossom became even
        //  Delta  - number of dual adjustments experience by vertices in the blossom as odd
        //           vertices before it became a root or even
        MaximumWeightMatching::intweight z0, t_root, t_odd, t_even, Delta;

        // Parent and children in the blossom tree
        Blossom* parent;
        std::list<std::pair<Blossom*, Edge>> subblossoms;

        // Status of the blossom in the search
        Label label;
        Edge backtrack_edge;

        // Flag needed for efficient backtracking
        bool visited;

        // Iterator to the blossom position in a shell's blossom list
        std::list<Blossom*>::iterator shell_blossoms_it;

        // For a non-even blossom maintains a list of it's nodes, the key for each node
        // is used to maintain a minimum slack of an edge connecting it to an even vertex
        SplitFindMin<NetworKit::node, MaximumWeightMatching::intweight, Edge, Blossom*>::List* list;

        bool is_trivial();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        std::list<NetworKit::node> node_list();
    };

    struct OldBlossom {
        // The number of all nodes in the old blossom
        NetworKit::count size;

        // Parent and child in heavy path decomposition
        OldBlossom* heavy_path_parent;
        OldBlossom* heavy_child;

        // Dual weight of the old blossom
        MaximumWeightMatching::intweight z;

        // Children of the shell
        std::list<OldBlossom*> children;

        // List of nodes solely in the shell
        std::list<NetworKit::node> nodes;

        // Blossoms contained in the shell
        std::list<Blossom*> shell_blossoms;

        // Sum of dual weights of the active shell and above at the time it became active
        MaximumWeightMatching::intweight current_shell_dual;

        // Status of the shell - wether it has been searched or dissolved
        bool dissolved, searched;

        // Index of the shell in the path
        NetworKit::index shell_index;

        bool is_heavy_path_root();
        void for_nodes(const std::function<void(NetworKit::node)>& handle);
        void for_blossoms(const std::function<void(OldBlossom*)>& handle);
    };

    // The working graph. It is the same as the provided graph when a perfect matching is sought.
    // Otherwise, it is the original graph reduced to an instance of MWPM.
    NetworKit::Graph workingGraph;
    std::vector<std::pair<NetworKit::node, NetworKit::node>> graph_edges;

    // Current edge weights and vertex dual variables
    std::vector<MaximumWeightMatching::intweight> current_w;
    std::vector<MaximumWeightMatching::intweight> current_y;

    // Trivial blossom for each node
    std::vector<Blossom*> trivial_blossom;

    // State of the current matching
    std::vector<NetworKit::node> matched_vertex;
    std::vector<NetworKit::edgeid> matched_edge;
    std::vector<bool> edge_in_matching;
    NetworKit::count edges_in_matching;

    void delete_blossom(Blossom* B, OldBlossom* S);
    void add_blossom(Blossom* B, OldBlossom* S);
    void expand_blossom(Blossom* B, OldBlossom* S);
    void swap_edge_in_matching(NetworKit::edgeid edge);
    void set_edge_in_matching(NetworKit::edgeid edge);
    void remove_edge_from_matching(NetworKit::edgeid edge);
    void check_edge_in_matching(NetworKit::edgeid edge);

    void create_trivial_blossoms(OldBlossom* T);
    OldBlossom* turn_current_blossoms_into_old(const std::list<Blossom*>& subblossoms);
    void heavy_path_decomposition(OldBlossom* T, MaximumWeightMatching::intweight outer_dual);
    void change_blossom_base(Blossom* B, NetworKit::node new_base, Edge edge);
    void swap_edges_on_even_path(Blossom* B, NetworKit::node u, NetworKit::node v);
    void swap_edges_on_even_path(Blossom* B, NetworKit::node out_vertex,
        std::list<Blossom*>&& out_blossoms);

    std::tuple<
        std::vector<NetworKit::node>,
        std::vector<MaximumWeightMatching::intweight>,
        OldBlossom*> scale(const std::vector<MaximumWeightMatching::intweight>& w);

    void match(OldBlossom* T);

    void path(OldBlossom* B, MaximumWeightMatching::intweight outer_dual);

    void shell_search(OldBlossom* B);
    void init_shell_search();

    NetworKit::count free_nodes_in_shell(OldBlossom* B);
    void enumerate_shells();
    void augmentPaths();
    MaximumWeightMatching::intweight current_slack(Edge e);

    void finish_shell_search();
    void delete_lists(Blossom* B);

    struct Event {
        enum Type { grow, blossom, dissolve, dissolveShell };
        Type type;
        union {
            struct GrowArgs { NetworKit::node v; Edge e; } grow;
            Edge uv;
            Blossom* B;
            OldBlossom* S;
        } args;

        static Event make_grow(NetworKit::node v, Edge e) {
            Event event;
            event.type = grow;
            event.args.grow.v = v; event.args.grow.e = e;
            return event;
        }

        static Event make_blossom(Edge uv) {
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

        static Event make_dissolveShell(OldBlossom* S) {
            Event event;
            event.type = dissolveShell;
            event.args.S = S;
            return event;
        }
    };

    // Outermost shell
    OldBlossom* old_root;

    // First shell in the current path
    OldBlossom* path_root;

    // First undissolved shell in the current path
    OldBlossom* highest_undissolved;

    // Shells in the current path and the number of free vertices it contains
    std::vector<std::pair<NetworKit::count, OldBlossom*>> shells;

    // Maps the vertices into the contracted graph during path augmentation
    std::vector<NetworKit::node> actual_to_contracted;

    // The blossom the vertex belongs to before augmentation
    std::vector<Blossom*> current_blossom;

    // Root of the path in which
    std::vector<OldBlossom*> vertex_path;

    // The shell the vertex belongs to before a path iteration
    std::vector<OldBlossom*> current_shell;

    // The status of the current search
    bool search_done;

    // Queue of search events
    ArrayPriorityQueue<Event> event_queue;

    // Time the old blossom has become the outer boundary of the shell in the search
    MaximumWeightMatching::intweight t_outer;

    // Time the old blossom has become the inner boundary in the search
    MaximumWeightMatching::intweight t_inner;

    // Time the whole graph has becoe the outer boundary in the search
    MaximumWeightMatching::intweight t_undissolvable;

    // Sum of the dual variables of blossoms above current path
    MaximumWeightMatching::intweight z_outer;

    // Sum of dual weights of the outer boundary and above at the time it became the boundary
    MaximumWeightMatching::intweight z_boundary;

    // Currently searched shell
    OldBlossom* active_shell;

    // Starting shell for the current shell
    OldBlossom* starting_shell;

    // The value of vertex's dual variable at the start of search
    std::vector<MaximumWeightMatching::intweight> y0;

    // The time the vertex was added to the current search
    std::vector<MaximumWeightMatching::intweight> t_shell;

    // Distributions to vertex's dual variable before it was searched
    std::vector<MaximumWeightMatching::intweight> Delta;

    // The starting shell of the search in which the vertex was searched
    std::vector<OldBlossom*> search_shell;

    // Used to maintain distributions by shells in the path
    FenwickTree shell_distribution;

    // Maintains nodes in even blossoms
    UnionFind<NetworKit::node, Blossom*> union_find;

    // Maintains the list for non-even blossoms
    SplitFindMin<NetworKit::node, MaximumWeightMatching::intweight, Edge, Blossom*> split_findmin;

    void grow(NetworKit::node v, Edge e);
    void schedule(NetworKit::node v);
    void dissolve(Blossom* B);
    void blossom(Edge edge);
    void dissolveShell(OldBlossom* S);

    struct BacktrackInfo {
        Blossom* blossom;
        Edge edge;
    };

    bool backtrack(Blossom* u, Blossom* v, Edge edge);
    bool backtrack_step(Blossom*& iter, std::vector<BacktrackInfo>& path);

    void create_new_blossom(Blossom* u, Blossom* v, Edge edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);

    std::list<Blossom*> blossoms_containing(NetworKit::node u, Blossom* until);
    void cut_path_at(std::vector<BacktrackInfo>& path, Blossom* cut, Blossom* cut_2);

    std::pair<std::list<std::pair<Blossom*, Edge>>, std::list<std::pair<Blossom*, Edge>>>
    split_subblossoms(std::list<std::pair<Blossom*, Edge>> subblossoms, Blossom* blossom);

    MaximumWeightMatching::intweight y(NetworKit::node v);
    MaximumWeightMatching::intweight z(Blossom* B);
    MaximumWeightMatching::intweight shell_z();
    MaximumWeightMatching::intweight slack(Edge e);
    bool is_exposed(Blossom* b);
    bool is_even(NetworKit::node v);
    bool is_odd(NetworKit::node v);
    Blossom* get_blossom(NetworKit::node v);
    void add_distribution(OldBlossom* S, MaximumWeightMatching::intweight distribution);
    MaximumWeightMatching::intweight distribution_so_far(int shell_index);
    bool matching_is_perfect();
    void clear_old_blossoms(OldBlossom* T);

    static std::string label_to_str(Label label);
    static Edge reverse(const Edge& edge);
    static void print_backtrack(
            Blossom* u, Blossom* v, Edge edge,
            std::vector<BacktrackInfo>& u_path, std::vector<BacktrackInfo>& v_path);
};

/**
 * @ingroup matching
 * The base class for the maximum cardinality matching algorithms.
 */
class MaximumCardinalityMatching :  public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the maximum matching algorithm.
     *
     * @param graph The input graph.
     */
    explicit MaximumCardinalityMatching(NetworKit::Graph &graph);

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
 * The class for the Micali-Vazirani maximum cardinality matching algorithm.
 */
class MicaliVaziraniMatching final :  public MaximumCardinalityMatching {
 public:
    explicit MicaliVaziraniMatching(NetworKit::Graph &graph);
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
