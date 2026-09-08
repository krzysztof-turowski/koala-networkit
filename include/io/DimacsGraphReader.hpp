/*
 * DimacsGraphReader.hpp
 *
 *  Created on: 21.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <unordered_map>

#include <networkit/io/GraphReader.hpp>

namespace Koala {

/**
 * @ingroup io
 * A reader for DIMACS graph format.
 * Full definition: http://prolland.free.fr/works/research/dsat/dimacs.html
 *
 */
class DimacsGraphReader final : public NetworKit::GraphReader {
 public:
    using McfResult = std::tuple<NetworKit::Graph,
        std::unordered_map<NetworKit::Edge, int64_t>,
        std::unordered_map<NetworKit::node, int64_t>,
        NetworKit::node, NetworKit::node>;
    using MinCostFlowResult = std::tuple<NetworKit::Graph,
        std::unordered_map<NetworKit::Edge, int64_t>,
        std::unordered_map<NetworKit::node, int64_t>>;

    DimacsGraphReader() = default;

    /**
     * Given the path of an input file, read the graph.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file
     */
    NetworKit::Graph read(std::string_view path) override;

    /**
     * Given the path of an input file, read the graph.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file, together with source and target nodes
     */
    std::tuple<NetworKit::Graph, NetworKit::node, NetworKit::node> read_all(
        const std::string &path);

    /**
     * Given the path of an input file,
     * read the graph together with minimum cost flow parameters.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file, together with source and target nodes
     */
    McfResult read_all_mcf(const std::string &path);

    /**
     * Given the path of an input file, read the graph for minimum cost flow.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file with edges,
     *              together with maps mapping edges to costs and nodes to supply/demand
     */
    MinCostFlowResult read_minimum_cost_flow(const std::string &path);
};

} /* namespace Koala */
