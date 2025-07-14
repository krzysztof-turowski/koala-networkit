#include <flow/MaximumFlow.hpp>

namespace Koala {

class PushRelabel final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;
    void run();
 private:
    int V;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> flow;
    std::unordered_map<NetworKit::node, int> height, nextedge, excess;
    std::queue<NetworKit::node> q;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> capacity;

    NetworKit::Edge reverse(const NetworKit::Edge &);

    void push(const NetworKit::node&, const NetworKit::Edge&);
    void relabel(const NetworKit::node&);
    void discharge(const NetworKit::node&);
    void initialize();
};

}  /* namespace Koala */
