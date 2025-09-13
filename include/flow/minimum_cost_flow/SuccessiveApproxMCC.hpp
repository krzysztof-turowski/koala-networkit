#pragma once

#include <set>
#include <unordered_set>
#include <list>

#include <flow/MinimumCostFlow.hpp>

namespace Koala {

class SuccessiveApproxMCC final : public MinimumCostFlow {
 public:   
    using MinimumCostFlow::MinimumCostFlow;

 private:
    void runImpl() override;
    void initFlow() override;
    void initCirculation() override;

    void initialize();
    void push(NetworKit::node, NetworKit::node, NetworKit::edgeid, bool reversed);
    void relabel(NetworKit::node const&);
    void refine();
    void wave();
    bool discharge(NetworKit::node const&);

    double cp(NetworKit::node, NetworKit::node, NetworKit::edgeid);
    int uf(NetworKit::edgeid, bool reversed);
    void force_flow(NetworKit::node, NetworKit::node, NetworKit::edgeid, int);
    node_map<double> potential;
    node_map<int> excess;
    node_map<int> pr_id;
    edgeid_map<int> lowerbound, upperbound;
    std::unordered_set<NetworKit::node> active;

    double epsi{0.};

    class DischargeList {
     public:
        virtual NetworKit::node getNext() { return 0; }
        virtual void moveToStart() {}
    };

    class ToposortList : public DischargeList {
     public:
        ToposortList(SuccessiveApproxMCC&);
        
        NetworKit::node getNext() override;
        void moveToStart() override;
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
