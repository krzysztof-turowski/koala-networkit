#ifndef FominKratschWoeginger_MDS_HPP_
#define FominKratschWoeginger_MDS_HPP_

#include <set>
#include <dominatingset/MinimumDominatingSet.hpp>

class FominKratschWoegingerMDS : public MinimumDominatingSet {
public:
    using MinimumDominatingSet::MinimumDominatingSet;
    FominKratschWoegingerMDS(const NetworKit::Graph &G);
    void run() override;
};

#endif /* EXHAUSTIVE_MDS_HPP_ */
