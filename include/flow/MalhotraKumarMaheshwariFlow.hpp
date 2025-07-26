#include <flow/MaximumFlow.hpp>

namespace Koala {

class MalhotraKumarMaheshwariFlow final : public MaximumFlow {
 public:
    using MaximumFlow::MaximumFlow;
    void run();
 private:
    int V;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> flow;
    NetworKit::Graph graph_stage;
    std::unordered_map<NetworKit::node, int> level;
    std::unordered_map<NetworKit::node, int> inPotential, outPotential;
    std::unordered_map<NetworKit::Edge, int, EdgeHash, EdgeEqual> capacity;

    NetworKit::Edge reverse(const NetworKit::Edge &);

    bool buildLevelGraph();
    void computePotential();
    void pushForward(NetworKit::node, NetworKit::edgeweight);
    void pushBackward(NetworKit::node, NetworKit::edgeweight);
    void deleteNode(NetworKit::node);
    void initialize();
};

}  /* namespace Koala */
