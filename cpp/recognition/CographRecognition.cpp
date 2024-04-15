#include <list>
#include "graph/GraphTools.hpp"
#include "recognition/CorneilStewartPerlCographRecognition.hpp"

namespace Koala {

    bool checkPath(const NetworKit::Graph &graph, NetworKit::node x, NetworKit::node y, NetworKit::node u, NetworKit::node v){
        return graph.hasEdge(y, u) && !graph.hasEdge(x, u) && !graph.hasEdge(x, v) && !graph.hasEdge(y, v);
    }
    void CographRecognition::check() const {
        assureFinished();
        for (const auto [x, y]: graph.edgeRange()) {
            for (const auto [u, v]: graph.edgeRange()) {
                if (x == u || x == v || y == u || y == v) {
                    continue;
                }
                assert(!checkPath(graph, x, y, u, v) && !checkPath(graph, x, y, v, u) && !checkPath(graph, y, x, u, v) && !checkPath(graph, y, x, v, u));
            }
        }
    }
}  // namespace Koala
