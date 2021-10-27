/*
 * DimacsBinaryGraphWriter.cpp
 *
 *  Created on: 27.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <algorithm>
#include <fstream>
#include <sstream>

#include <networkit/auxiliary/Enforce.hpp>

#include <io/DimacsBinaryGraphWriter.hpp>

namespace Koala {

void DimacsBinaryGraphWriter::write(const NetworKit::Graph &G, const std::string &path) {
    std::ofstream graphFile(path);
    Aux::enforceOpened(graphFile);
    
    std::stringstream preambleFile;
    preambleFile << "p edge " << G.numberOfNodes() << ' ' << G.numberOfEdges() << '\n';
    std::string preamble(preambleFile.str());
    graphFile << preamble.size() << std::endl << preamble;

    std::string data((G.numberOfNodes() >> 3) + 1, 0);
    const unsigned char MASKS[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
    G.forNodes([&](NetworKit::node v) {
        std::fill_n(data.data(), (v >> 3) + 1, 0);
        G.forNeighborsOf(v, [&](NetworKit::node u) {
            if (u <= v) {
                data[u >> 3] |= MASKS[u & 7];
            }
        });
        graphFile.write(data.data(), (v >> 3) + 1);
    });
}

} /* namespace Koala */
