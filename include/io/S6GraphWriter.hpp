/*
 * G6GraphWriter.hpp
 *
 *  Created on: 20.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <string>

#include <networkit/io/GraphWriter.hpp>

namespace Koala {

/**
 * @ingroup io
 * A writer for sparse6 graph format. Each line contains a single graph.
 * Full definition: https://users.cecs.anu.edu.au/~bdm/data/formats.txt
 *
 */
class S6GraphWriter final : public NetworKit::GraphWriter {
 public:
    S6GraphWriter() = default;

    /**
     * Given a graph and a file path, write the graph to the file.
     *
     * @param[in]  G     input graph
     * @param[in]  path  output file path
     */
    void write(const NetworKit::Graph &G, const std::string &path) override;

    /**
     * Given a graph, find its graph6 representation.
     *
     * @param[in]  G     input graph
     * @param[out]  output string
     */
    std::string writeline(const NetworKit::Graph &G);
};

} /* namespace Koala */
