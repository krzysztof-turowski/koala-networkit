#include <flow/MaximumFlow.hpp>

namespace Koala {

class BoykovKolmogorovFlow final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;
    void run();
 private:
    int V;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> flow;
    NetworKit::node spath, tpath;
    std::unordered_map<NetworKit::node, NetworKit::node> parent;
    std::unordered_map<NetworKit::node, int> tree;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> capacity;
    std::queue<NetworKit::node> active;
    std::queue<NetworKit::node> orphan;

    NetworKit::Edge reverse(const NetworKit::Edge &);

    int tree_capacity(NetworKit::node, NetworKit::node);
    void initialize();
    bool grow();
    int augment();
    void adopt();
    bool origin(NetworKit::node);
};

}  /* namespace Koala */
