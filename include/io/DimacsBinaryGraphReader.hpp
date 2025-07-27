/*
 * DimacsBinaryGraphReader.hpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <string>

#include <networkit/io/GraphReader.hpp>

namespace Koala {

/**
 * @ingroup io
 * A reader for DIMACS binary graph format.
 * Full definition: https://mat.tepper.cmu.edu/COLOR/format/README.binformat
 *
 */
class DimacsBinaryGraphReader final : public NetworKit::GraphReader {
 public:
    DimacsBinaryGraphReader() = default;

    /**
     * Given the path of an input file, read the graph.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file
     */
    NetworKit::Graph read(std::string_view path) override;
};

}  /* namespace Koala */
