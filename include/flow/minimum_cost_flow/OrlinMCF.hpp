#pragma once

#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Attributes.hpp>

namespace Koala {

class OrlinMCF : public MinimumCostFlow {
 private:
    const double ALPHA = 0.9;
    int imbalanced = 0;
    int delta = 0;


    void runImpl() override;
    void initialize();
    void setupValues();
    void makeConnected();
    void contractNodes(NetworKit::Graph&, NetworKit::node, NetworKit::node);
    void deltaScalingPhase(NetworKit::Graph&);
    void contractionPhase(NetworKit::Graph&);

    NetworKit::Graph::NodeIntAttribute base(NetworKit::Graph &);
    NetworKit::Graph::NodeIntAttribute excess(NetworKit::Graph &);
    NetworKit::Graph::NodeIntAttribute potential(NetworKit::Graph &g);
    NetworKit::Graph::EdgeIntAttribute flow(NetworKit::Graph &g);

};

} /* namespace Koala */

 