#pragma once

#include <vector>

#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * Electrical network solver used by ElectricalFlow.
 *
 * Given vertex demands and edge resistances, compute() obtains Laplacian potentials and the
 * induced electrical flow.
 */
class ElectricalNetwork {
 public:
    /**
     * Construct an electrical network for a fixed graph and demand vector.
     *
     * @param G Input graph.
     * @param demand Net demand at each vertex.
     */
    ElectricalNetwork(const NetworKit::Graph &G, const std::vector<double> &demand);

    /**
     * Compute potentials and induced flows for the supplied edge resistances.
     *
     * @param resistance Resistance matrix indexed by edge endpoints.
     */
    void compute(const std::vector<std::vector<double>> &resistance);

    std::vector<std::vector<double>> flow;
    std::vector<double> potentials;

 private:
    const NetworKit::Graph &graph;
    const std::vector<double> &demand;
};
}  // namespace Koala
