/*
 * DimacsBinaryGraphWriter.hpp
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
 * A writer for DIMACS binary graph format.
 * Full definition: https://mat.tepper.cmu.edu/COLOR/format/README.binformat
 *
 */
class DimacsBinaryGraphWriter final : public NetworKit::GraphWriter {

public:
    DimacsBinaryGraphWriter() = default;

    /**
     * Given a graph and a file path, write the graph to the file in DIMACS binary format.
     *
     * @param[in]  G     input graph
     * @param[in]  path  output file path
     */
    void write(const NetworKit::Graph &G, const std::string &path) override;
};

} /* namespace Koala */
