/*
 * BranchAndReduceSetCover.hpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#pragma once

#include <set>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup set_cover
 * The class for exact determination of minimum set cover
 * via branch and reduce algorithm.
 *
 */
class BranchAndReduceSetCover : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the branch and reduce minimum set cover algorithm.
     *
     * @param family The input family of subsets.
     * @param occurences The number of occurences of each element.
     */
    BranchAndReduceSetCover(
        std::vector<std::set<NetworKit::node>> &family,
        std::vector<std::set<NetworKit::index>> &occurences);

    /**
     * Execute the minimum set cover procedure.
     */
    void run();

    /**
     * Return the minimum set cover found by the algorithm.
     *
     * @return Minimum set cover vector.
     */
    std::vector<bool> getSetCover() const;

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

 protected:
    std::vector<std::set<NetworKit::node>> &family;
    std::vector<std::set<NetworKit::index>> &occurences;
    std::vector<bool> set_cover;

    bool reduce();
    NetworKit::index find_unique_occurence_set();
    virtual bool reduce_matching() = 0;
    std::vector<bool> forced_set_cover(NetworKit::index index);
    std::vector<bool> discarded_set_cover(NetworKit::index index);
};

class GrandoniSetCover : public BranchAndReduceSetCover {
 public:
    using BranchAndReduceSetCover::BranchAndReduceSetCover;

 protected:
    virtual bool reduce_matching();
};

class FominGrandoniKratschSetCover : public BranchAndReduceSetCover {
 public:
    using BranchAndReduceSetCover::BranchAndReduceSetCover;

 protected:
    virtual bool reduce_matching();
};

class RooijBodlaenderSetCover : public FominGrandoniKratschSetCover {
 public:
    using FominGrandoniKratschSetCover::FominGrandoniKratschSetCover;

 private:
    bool reduce();
    NetworKit::index find_counting_rule_reduction_set();
    NetworKit::index find_cardinality_frequency_set();
};

}  /* namespace Koala */
