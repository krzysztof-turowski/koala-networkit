#ifndef ROOIJ_BODLAENDER_MDS_HPP_
#define ROOIJ_BODLAENDER_MDS_HPP_

#include <set>
#include <dominatingset/MinimumDominatingSet.hpp>

class RooijBodlaenderMDS : public MinimumDominatingSet {
public:
    using MinimumDominatingSet::MinimumDominatingSet;
    RooijBodlaenderMDS(const NetworKit::Graph &G);
    void run() override;
};

#endif /* ROOIJ_BODLAENDER_MDS_HPP_ */
