#pragma once

#include <set>
#include <unordered_set>
#include <list>

#include <flow/MinimumCostFlow.hpp>

namespace Koala {

class SuccessiveApproxMCC : public MinimumCostFlow {
 public:   
    using MinimumCostFlow::MinimumCostFlow;
    SuccessiveApproxMCC::SuccessiveApproxMCC(NetworKit::Graph&, 
        edge_map<MCFEdgeParams>&, node_map<int>&);
    SuccessiveApproxMCC::SuccessiveApproxMCC(NetworKit::Graph&, 
        edge_map<MCFEdgeParams>&, NetworKit::node, NetworKit::node, int);

 private:
    void runImpl();
    void initialize();
    bool push_relabel(NetworKit::node const&);
    void push(NetworKit::node const&, NetworKit::node const&);
    void relabel(NetworKit::node const&);
    void refine();
    void wave();
    bool discharge(NetworKit::node const&);

    double cp(NetworKit::node const&, NetworKit::node const&);
    int uf(NetworKit::node const&, NetworKit::node const&);
    void force_flow(NetworKit::node const&, NetworKit::node const&, int);
    node_map<double> price;
    node_map<int> excess;
    node_map<int> pr_id;
    std::unordered_set<NetworKit::node> active;

    double epsi{0.};

    class DischargeList {
     public:
        virtual NetworKit::node getNext();
        virtual void moveToStart();
    };

    class ToposortList : public DischargeList {
        ToposortList(SuccessiveApproxMCC&);
        
        NetworKit::node getNext();
        void moveToStart();
     private:
        SuccessiveApproxMCC &approx;
        std::list<NetworKit::node> nodes;
        node_map<bool> vis; 
        void dfs(NetworKit::node);

        std::list<NetworKit::node>::iterator it1;
        std::list<NetworKit::node>::iterator it2;
    };

};

} /* namespace Koala */
