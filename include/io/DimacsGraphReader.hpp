/*
 * DimacsGraphReader.hpp
 *
 *  Created on: 21.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <string>

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
    DimacsGraphReader() = default;

    /**
     * Given the path of an input file, read the graph.
     *
     * @param[in]  path  input file path
     * @param[out]  the graph read from file
     */
    NetworKit::Graph read(const std::string &path) override;
};

} /* namespace Koala */
