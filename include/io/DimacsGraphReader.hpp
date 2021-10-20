/*
 * DimacsGraphReader.hpp
 *
 *  Created on: 21.10.2021
 *      Author: Krzysztof Turowski
 */

#ifndef KOALA_IO_DIMACS_GRAPH_READER_HPP_
#define KOALA_IO_DIMACS_GRAPH_READER_HPP_

#include <string>

#include <networkit/io/GraphReader.hpp>

namespace Koala {

/**
 * @ingroup io
 * A reader for DIMACS graph format. Each line contains a single graph.
 * Full definition: http://archive.dimacs.rutgers.edu/pub/challenge/graph/doc/ccformat.tex
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
#endif // KOALA_IO_DIMACS_GRAPH_READER_HPP_
