#include <iostream>
#include <fstream>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <vector>
#include <chrono>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "../bakers_technique/bakers_technique.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kouter
#include <boost/test/unit_test.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace boost;

struct file_reader{
    std::ifstream file;
    std::string edges;

    file_reader(std::string fn):
        file("test_graphs/" + fn){}

    bool next_graph(Graph& g) {
        if(!file.is_open()) {
            return false;
        }

        int n, m;

        if (!(file >> n >> m)) {
            return false;
        }

        edges = std::to_string(n) + " " + std::to_string(m) + "\n";

        int a, b;
        while (m--) {
            file >> a >> b;
            add_edge(a, b, g);
            edges += std::to_string(a) + " " + std::to_string(b) + "  ";
        }

        edges += "\n";

        return true;
    }

    bool next_graph2(Graph& g) {
        if(!file.is_open()) {
            return false;
        }
        edges.clear();
        std::string s;
        std::set <std::pair<int, int> > es;
        while (getline(file, s)) {
            if (s.empty()) {
                break;
            }
            else {
                std::istringstream tmp(s);
                int a, b;
                char nic;
                tmp >> a >> nic >> b;
                a--; b--;
                std::pair<int, int> e(std::min(a,b), std::max(a,b));
                if (a != b && es.find(e) == es.end()) {
                    add_edge(a, b, g);
                    es.insert(e);
                    edges += std::to_string(a) + " " + std::to_string(b) + "  ";
                }
            }
        }
        edges = std::to_string(num_vertices(g)) + " " + std::to_string(num_edges(g)) + "\n" + edges;
        return true;
    }

    std::string get_last_edges(){
        return edges;
    }
};

void make_graph(Graph& g, int m, std::string edges) {
    std::stringstream stream(edges);
    int a, b;
    while (m--) {
        stream >> a >> b;
        add_edge(a, b, g);
    }
}

void get_embedding(Graph& g, PlanarEmbedding& embedding, std::vector<int>& outer_face) {
    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);


    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )) {}
//        std::cout << "Input graph is planar" << std::endl;
    else {
        std::cout << "Input graph is not planar" << std::endl;
        return;
    }

    std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    face_getter<Edge> my_vis(&faces, vertices_in_face);
    level_face_traversal<Graph>(embedding, my_vis);

    for (int v : vertices_in_face[0]) {
        outer_face.push_back(v);
    }
}

int independent_set_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mx = 0;

    int last = n;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (s[e.m_source] && s[e.m_target]) {
                res = false;
                break;
            }
        }

        if (res) {
            int ones = s.count();
            mx = std::max(mx, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }

//        int f = s.find_first();
//        if (f < last) {
//            last = f;
//            std::cout << f << std::endl;
//        }
    }

    return mx;
}

int vertex_cover_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mn = INT16_MAX;

    int last = n;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (!s[e.m_source] && !s[e.m_target]) {
                res = false;
                break;
            }
        }

        if (res) {
            int ones = s.count();
            mn = std::min(mn, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }
    }

    return mn;
}

int dominating_set_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mn = INT16_MAX;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        boost::dynamic_bitset<> dom(n, 0);
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (s[e.m_source] || s[e.m_target]) {
                dom[e.m_source] = true;
                dom[e.m_target] = true;
            }
        }

        if (dom.count() == n) {
            int ones = s.count();
            mn = std::min(mn, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }
    }

    return mn;
}

BOOST_AUTO_TEST_SUITE(corectness)
    BOOST_AUTO_TEST_CASE(przyklad) {
        Graph g;
        make_graph(g, 158, "0 1  0 4  1 2  2 5  2 106  2 110  2 4  2 3  3 7  3 22  3 21  3 11  3 8  3 4  4 5  5 6  6 7  7 8  7 105  8 9  8 103  9 10  10 11  10 13  10 14  10 19  10 20  10 97  10 25  10 23  10 103  11 12  13 17  13 14  14 15  14 16  14 17  16 17  17 18  20 101  20 21  21 22  21 25  21 94  21 93  21 91  21 95  21 34  21 37  21 36  21 96  21 100  22 23  23 24  25 26  25 27  25 28  25 89  25 47  25 46  25 41  25 34  25 90  25 93  28 29  28 86  28 30  28 84  28 31  28 87  28 88  29 30  30 83  30 85  30 31  31 82  31 32  32 46  32 33  33 38  33 39  33 40  33 52  33 81  33 34  34 35  34 36  36 37  38 39  40 56  40 55  40 80  40 41  41 43  41 61  41 42  43 65  43 63  43 79  43 62  43 44  44 45  44 64  45 69  45 64  45 78  45 46  46 47  46 51  46 52  46 77  46 56  46 70  47 48  47 49  47 50  52 53  53 54  53 55  53 56  56 57  57 58  57 59  59 60  59 73  59 69  60 61  60 62  62 63  62 64  64 67  64 68  64 69  64 71  64 72  64 76  64 65  65 66  68 69  69 71  69 70  72 73  73 74  73 75  83 84  88 89  90 91  91 92  96 97  96 98  96 99  97 98  100 101  101 102  103 104  106 107  107 108  107 109");
        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);
        int result = baker<independent_set>(g, embedding, outer_face);
    }

    BOOST_AUTO_TEST_CASE(is_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(is_five_vertices) {
        file_reader f("5vertices");
//        std::string s = get_current_dir_name();

        int i = 1;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(is_six_vertices) {
        file_reader f("6vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(is_seven_vertices) {
        file_reader f("7vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(is_eight_vertices) {
        file_reader f("8vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(is_nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_five_vertices) {
        file_reader f("5vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_six_vertices) {
        file_reader f("6vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_seven_vertices) {
        file_reader f("7vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_eight_vertices) {
        file_reader f("8vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_five_vertices) {
        file_reader f("5vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_six_vertices) {
        file_reader f("6vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_seven_vertices) {
        file_reader f("7vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_eight_vertices) {
        file_reader f("8vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
            if (result != expected) {
                break;
            }
        }
    }

    BOOST_AUTO_TEST_CASE(tr_simple)
    {
        Graph g;
        make_graph(g, 3, "0 2  1 3  2 3");

        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);

        bodlaender_dominating_set_ecc(g, embedding, outer_face);
    }

    BOOST_AUTO_TEST_CASE(tr_vc_seven_vertices) {
        file_reader f("7vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_vc_eight_vertices) {
        file_reader f("8vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_vc_nine_vertices) {
        file_reader f("9vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_seven_vertices) {
        file_reader f("7vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_eight_vertices) {
        file_reader f("8vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_nine_vertices) {
        file_reader f("9vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_lcc_seven_vertices) {
        file_reader f("7vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_dominating_set_lcc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_lcc_eight_vertices) {
        file_reader f("8vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_dominating_set_lcc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_ecc_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_dominating_set_ecc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_ecc_six_vertices) {
        file_reader f("6vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_dominating_set_ecc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }
BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(parsers)
    BOOST_AUTO_TEST_CASE(graph_sorter) {
        file_reader f("unprocessed_graphs/same_low_n");

        int z = 10000;
        bool res = true;
        std::vector<std::ofstream> files(15);
        for (int i = 1; i < 15; i++) {
            files[i].open("test_graphs/performance_test_graphs/" + std::to_string(i) + "-outer");
            if (!files[i].is_open()) {
                return;
            }
        }
        std::vector<int> graphs(10);
        while (z--) {
            Graph g;
            res = f.next_graph2(g);
            if (num_vertices(g) == 0) {
                continue;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            std::vector<int> vertex_level(num_vertices(g));
            std::vector< std::vector<Edge> > aaaa;
            int level = name_levels(embedding, outer_face, vertex_level, aaaa);
            if (level < 10)
                files[level] << f.get_last_edges() << "\n";
        }
    }

    BOOST_AUTO_TEST_CASE(graph_sorter2) {
        file_reader f("unprocessed_graphs/different_n");
        std::ofstream file_small("test_graphs/performance_test_graphs/small_graphs");
        std::ofstream file_big("test_graphs/performance_test_graphs/big_graphs");
        int z = 100;
        bool res = true;
        std::vector< std::vector<std::string> > graphs(10000);
        std::vector<int> num;
        while (z--) {
            Graph g;
            res = f.next_graph2(g);
            if (num_vertices(g) == 0) {
                continue;
            }
            num.push_back(num_vertices(g));
            graphs[num_vertices(g)].push_back(f.get_last_edges());
        }

        std::sort(num.begin(), num.end());
        int max_size = num[num.size()/2];

        for (int n = 0; n < max_size; n++) {
            auto& gs = graphs[n];
            for (auto& g : gs) {
                file_small << g << "\n";
            }
        }
        for (int n = max_size; n < 1000; n++) {
            auto& gs = graphs[n];
            for (auto& g : gs) {
                file_big << g << "\n";
            }
        }
    }
BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(performance)
    BOOST_AUTO_TEST_CASE(kouter_is_performance) {
        bool res = true;
        std::vector<double> results_baker(15);
        std::vector<double> results_bodlaender(15);
        std::vector<int> num_of_graphs(15);
        std::cout << "k: " << std::flush;
        for (int i = 3; i < 10; i++) {
            file_reader f("performance_test_graphs/" + std::to_string(i) + "-outer");
            int z=10;
            std::cout << i << " g: " << std::flush;
            while (z--) {
                std::cout << z << std::flush;
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);

                auto start = std::chrono::steady_clock::now();
                baker<independent_set>(g, embedding, outer_face);
                auto stop = std::chrono::steady_clock::now();
                results_baker[i] += std::chrono::duration<double, std::milli>(stop - start).count();

                start = std::chrono::steady_clock::now();
                bodlaender_independent_set(g, embedding, outer_face);
                stop = std::chrono::steady_clock::now();
                results_bodlaender[i] += std::chrono::duration<double, std::milli>(stop - start).count();

                num_of_graphs[i]++;

                std::cout << "\b" << std::flush;
            }

            std::cout << "\b\b\b\b\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        std::ofstream file_baker("results/k-outer_baker_performance");
        file_baker << "Baker\n";
        std::ofstream file_bodlaender("results/k-outer_bodlaender_performance");
        file_bodlaender << "Bodlaender\n";
        for (int i = 3; i < 10; i++) {
            if (num_of_graphs[i] > 0) {
                file_baker << i << " " << results_baker[i] / num_of_graphs[i] << "\n";
                file_bodlaender << i << " " << results_bodlaender[i] / num_of_graphs[i] << "\n";
            }
        }
        file_baker.close();
        file_bodlaender.close();

        system("python3 ../plot/plot.py k-outer_baker_performance k-outer_bodlaender_performance is_k.pdf \"zewnętrzna planarność\"");
    }

    BOOST_AUTO_TEST_CASE(n_is_performance) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        file_reader f("performance_test_graphs/small_graphs");
        std::ofstream file_baker("results/n_baker_performance");
        file_baker << "Baker\n";
        std::ofstream file_bodlaender("results/n_bodlaender_performance");
        file_bodlaender << "Bodlaender\n";

        std::cout << "g: " << std::flush;

        int z = 250;
        while (z--) {
            std::cout << z << std::flush;
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);

            auto start = std::chrono::steady_clock::now();
            baker<independent_set>(g, embedding, outer_face);
            auto stop = std::chrono::steady_clock::now();
            file_baker << num_vertices(g) << " " << std::chrono::duration<double, std::milli>(stop - start).count() << "\n";

            start = std::chrono::steady_clock::now();
            bodlaender_independent_set(g, embedding, outer_face);
            stop = std::chrono::steady_clock::now();
            file_bodlaender << num_vertices(g) << " " << std::chrono::duration<double, std::milli>(stop - start).count() << "\n";

            for (int _ = 0; _ < std::to_string(z).size(); _++)
                std::cout << "\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        file_baker.close();
        file_bodlaender.close();

        system("python3 ../plot/plot.py n_baker_performance n_bodlaender_performance is_n.pdf \"liczba wierzchołków\"");
    }

    BOOST_AUTO_TEST_CASE(ptas_is_performance) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        std::ofstream file_baker("results/ptas_baker_performance");
        file_baker << "Baker\n";
        std::ofstream file_bodlaender("results/ptas_bodlaender_performance");
        file_bodlaender << "Bodlaender\n";
        std::cout << "k: " << std::flush;
        for (int k = 1; k < 10; k++) {
            std::cout << k << " g: " << std::flush;

            file_reader f("performance_test_graphs/big_graphs");
            int z = 10;
            double time_baker = 0, time_bodlaender = 0;
            while (z--) {
                std::cout << z << std::flush;
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);
                auto start = std::chrono::steady_clock::now();
                bakers_technique(g, embedding, outer_face, k, Baker, is);
                auto stop = std::chrono::steady_clock::now();
                time_baker += std::chrono::duration<double, std::milli>(stop - start).count();

                start = std::chrono::steady_clock::now();
                bakers_technique(g, embedding, outer_face, k, Bodlaender, is);
                stop = std::chrono::steady_clock::now();
                time_bodlaender += std::chrono::duration<double, std::milli>(stop - start).count();

                std::cout << "\b" << std::flush;
            }
            file_baker << k << " " << time_baker << "\n";
            file_bodlaender << k << " " << time_bodlaender << "\n";

            std::cout << "\b\b\b\b\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        file_baker.close();
        file_bodlaender.close();

        system("python3 ../plot/plot.py ptas_baker_performance ptas_bodlaender_performance is_ptas.pdf \"k\"");
    }

    BOOST_AUTO_TEST_CASE(n_ds_performance) {
        bool res = true;
        std::vector< std::pair<double, int> > results_baker(400);
        std::vector< std::pair<double, int> > results_bod_ecc(400);
        std::vector< std::pair<double, int> > results_bod_lcc(400);
        file_reader f("performance_test_graphs/small_graphs");
        std::ofstream file_baker("results/n_ds_baker_performance");
        file_baker << "Baker\n";
        std::ofstream file_ecc("results/n_ds_ecc_performance");
        file_ecc << "Bodlaender ECC\n";
        std::ofstream file_lcc("results/n_ds_lcc_performance");
        file_lcc << "Bodlaender LCC\n";

        std::cout << "g: " << std::flush;

        int z = 120;
        while (z--) {
            std::cout << z << std::flush;
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);

            auto start = std::chrono::steady_clock::now();
            baker<dominating_set>(g, embedding, outer_face);
            auto stop = std::chrono::steady_clock::now();
            results_baker[num_vertices(g)].first += std::chrono::duration<double, std::milli>(stop - start).count();
            results_baker[num_vertices(g)].second++;

            if (z > 100) {
                start = std::chrono::steady_clock::now();
                bodlaender_dominating_set_ecc(g, embedding, outer_face);
                stop = std::chrono::steady_clock::now();
                results_bod_ecc[num_vertices(g)].first += std::chrono::duration<double, std::milli>(stop - start).count();
                results_bod_ecc[num_vertices(g)].second++;

                start = std::chrono::steady_clock::now();
                bodlaender_dominating_set_lcc(g, embedding, outer_face);
                stop = std::chrono::steady_clock::now();
                results_bod_lcc[num_vertices(g)].first += std::chrono::duration<double, std::milli>(stop - start).count();
                results_bod_lcc[num_vertices(g)].second++;

            }

            for (int _ = 0; _ < std::to_string(z).size(); _++)
                std::cout << "\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        for (int i = 0; i < 400; i++) {
            if (results_baker[i].second > 0) {
                file_baker << i << " " << results_baker[i].first / results_baker[i].second << "\n";
            }
        }

        for (int i = 0; i < 400; i++) {
            if (results_bod_ecc[i].second > 0) {
                file_ecc << i << " " << results_bod_ecc[i].first / results_bod_ecc[i].second << "\n";
            }
        }

        for (int i = 0; i < 400; i++) {
            if (results_bod_lcc[i].second > 0) {
                file_lcc << i << " " << results_bod_lcc[i].first / results_bod_lcc[i].second << "\n";
            }
        }


        file_baker.close();
        file_ecc.close();
        file_lcc.close();

        system("python3 ../plot/plot.py n_ds_baker_performance n_ds_ecc_performance n_ds_lcc_performance ds_n.pdf \"liczba wierzchołków\"");
    }

    BOOST_AUTO_TEST_CASE(kouter_baker_performance) {
        bool res = true;
        std::vector<double> results_baker_is(15);
        std::vector<double> results_baker_ds(15);
        std::vector<int> num_of_graphs(15);
        std::cout << "k: " << std::flush;
        for (int i = 3; i < 10; i++) {
            file_reader f("performance_test_graphs/" + std::to_string(i) + "-outer");
            int z=5;
            std::cout << i << " g: " << std::flush;
            while (z--) {
                std::cout << z << std::flush;
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);

                auto start = std::chrono::steady_clock::now();
                baker<independent_set>(g, embedding, outer_face);
                auto stop = std::chrono::steady_clock::now();
                results_baker_is[i] += std::chrono::duration<double, std::milli>(stop - start).count();

                if (i < 7) {
                    start = std::chrono::steady_clock::now();
                    baker<dominating_set>(g, embedding, outer_face);
                    stop = std::chrono::steady_clock::now();
                    results_baker_ds[i] += std::chrono::duration<double, std::milli>(stop - start).count();
                }

                num_of_graphs[i]++;

                std::cout << "\b" << std::flush;
            }

            std::cout << "\b\b\b\b\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        std::ofstream file_baker_is("results/k-outer_baker_performance");
        file_baker_is << "zbiór niezależny\n";
        std::ofstream file_baker_ds("results/k-outer_ds_baker_performance");
        file_baker_ds << "zbiór dominujący\n";
        for (int i = 3; i < 10; i++) {
            if (num_of_graphs[i] > 0) {
                file_baker_is << i << " " << results_baker_is[i] / num_of_graphs[i] << "\n";
                file_baker_ds << i << " " << results_baker_ds[i] / num_of_graphs[i] << "\n";
            }
        }
        file_baker_is.close();
        file_baker_ds.close();

        system("python3 ../plot/plot.py k-outer_baker_performance k-outer_ds_baker_performance is_ds_k.pdf \"zewnętrzna planarność\"");
    }

    BOOST_AUTO_TEST_CASE(n_baker_performance) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        file_reader f("performance_test_graphs/small_graphs");
        std::ofstream file_baker_is("results/n_baker_performance");
        file_baker_is << "zbiór niezależny\n";
        std::ofstream file_baker_ds("results/n_ds_baker_performance");
        file_baker_ds << "zbiór dominujący\n";

        std::cout << "g: " << std::flush;

        int z = 200;
        while (z--) {
            std::cout << z << std::flush;
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);

            auto start = std::chrono::steady_clock::now();
            baker<independent_set>(g, embedding, outer_face);
            auto stop = std::chrono::steady_clock::now();
            file_baker_is << num_vertices(g) << " " << std::chrono::duration<double, std::milli>(stop - start).count() << "\n";

            start = std::chrono::steady_clock::now();
            baker<dominating_set>(g, embedding, outer_face);
            stop = std::chrono::steady_clock::now();
            file_baker_ds << num_vertices(g) << " " << std::chrono::duration<double, std::milli>(stop - start).count() << "\n";

            for (int _ = 0; _ < std::to_string(z).size(); _++)
                std::cout << "\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        file_baker_is.close();
        file_baker_ds.close();

        system("python3 ../plot/plot.py n_baker_performance n_ds_baker_performance is_ds_n.pdf \"liczba wierzchołków\"");
    }

    BOOST_AUTO_TEST_CASE(ptas_is_results) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        std::ofstream file("results/ptas_is_baker_results");

        std::cout << "k: " << std::flush;

        for (int k = 1; k < 10; k++) {
            std::cout << k << " g: " << std::flush;

            file_reader f("performance_test_graphs/big_graphs");
            int z = 10;
            double result = 0;
            while (z--) {
                std::cout << z << std::flush;
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);
                result += bakers_technique(g, embedding, outer_face, k, Baker, is);

                std::cout << "\b" << std::flush;
            }
            file << k << " " << result << "\n";

            std::cout << "\b\b\b\b\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        file.close();

        system("python3 ../plot/plot2.py results/ptas_is_baker_results max is_ptas_res_pdf");
    }

    BOOST_AUTO_TEST_CASE(ptas_vc_results) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        std::ofstream file("results/ptas_vc_baker_results");

        std::cout << "k: " << std::flush;

        for (int k = 1; k < 10; k++) {
            std::cout << k << " g: " << std::flush;

            file_reader f("performance_test_graphs/big_graphs");
            int z = 10;
            double result = 0;
            while (z--) {
                std::cout << z << std::flush;
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);
                result += bakers_technique(g, embedding, outer_face, k, Baker, vc);

                std::cout << "\b" << std::flush;
            }
            file << k << " " << result << "\n";

            std::cout << "\b\b\b\b\b" << std::flush;
        }

        std::cout << "\b\b\b" << std::flush;

        file.close();

        system("python3 ../plot/plot2.py results/ptas_vc_baker_results min vc_ptas_res.pdf");
    }

BOOST_AUTO_TEST_SUITE_END()

