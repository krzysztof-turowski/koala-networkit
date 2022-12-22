/*
 * DimacsGraphWriter.hpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <string>

#include <networkit/io/GraphWriter.hpp>

namespace Koala {

/**
 * @ingroup io
 * A writer for DIMACS graph format.
 * Full definition: http://prolland.free.fr/works/research/dsat/dimacs.html
 *
 */
class DimacsGraphWriter final : public NetworKit::GraphWriter {

public:
    DimacsGraphWriter() = default;

    /**
     * Given a graph and a file path, write the graph to the file in DIMACS format.
     *
     * @param[in]  G     input graph
     * @param[in]  path  output file path
     */
    void write(const NetworKit::Graph &G, const std::string &path) override;
};

} /* namespace Koala */
