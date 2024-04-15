#include <list>
#include "graph/GraphTools.hpp"
#include "recognition/CorneilStewartPerlCographRecognition.h"

namespace Koala {
    CorneilStewartPerlCographRecognition::CorneilStewartPerlCographRecognition(NetworKit::Graph &graph) : graph(graph), is_cograph(State::UNKNOWN) {

    }

    bool CorneilStewartPerlCographRecognition::isCograph() const {
        assureFinished();
        return is_cograph == State::COGRAPH;
    }

    CorneilStewartPerlCographRecognition::State CorneilStewartPerlCographRecognition::getState() const {
        assureFinished();
        return is_cograph;
    }

    CorneilStewartPerlCographRecognition::State recognition(NetworKit::Graph &graph) {
        return CorneilStewartPerlCographRecognition::Cograph_Recognition(graph);
    }

    void CorneilStewartPerlCographRecognition::run() {
        hasRun = true;
        if (is_cograph != State::UNKNOWN) {
            return;
        }
        is_cograph = recognition(graph);
        if (is_cograph != State::UNKNOWN) {
            return;
        }
        is_cograph = State::COGRAPH;
    }
    bool checkPath(const NetworKit::Graph &graph, NetworKit::node x, NetworKit::node y, NetworKit::node u, NetworKit::node v){
        return graph.hasEdge(y, u) && !graph.hasEdge(x, u) && !graph.hasEdge(x, v) && !graph.hasEdge(y, v);
    }
    void CorneilStewartPerlCographRecognition::check() const {
        assureFinished();
        for (auto e1: graph.edgeRange()) {
            for (auto e2: graph.edgeRange()) {
                auto [x, y] = e1;
                auto [u, v] = e2;
                if (x == u || x == v || y == u || y == v) {
                    continue;
                }
                assert(!checkPath(graph, x, y, u, v) && !checkPath(graph, x, y, v, u) && !checkPath(graph, y, x, u, v) && !checkPath(graph, y, x, v, u));
            }
        }
    }
}  // namespace Koala
