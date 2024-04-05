#include <list>
#include "graph/GraphTools.hpp"
#include "recognition/CorneilStewartPerlCographRecognition.h"
namespace Koala {
    CorneilStewartPerlCographRecognition::CorneilStewartPerlCographRecognition(NetworKit::Graph &graph):  graph(graph), is_cograph(State::UNKNOWN) {

    }
    bool CorneilStewartPerlCographRecognition::isCograph() const{
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

    void CorneilStewartPerlCographRecognition::check() const {
        assureFinished();
        bool is_cograph = true;
        for (auto e1 : graph.edgeRange()) {
            for (auto e2 : graph.edgeRange()) {
                auto x = e1.u;
                auto y = e1.v;
                auto u = e2.u;
                auto v = e2.v;
                if (x == u || x == v || y == u || y == v) {
                    continue;
                }
                if (graph.hasEdge(y, u) && !graph.hasEdge(x, u) && !graph.hasEdge(x, v) && !graph.hasEdge(y, v)) {
                    is_cograph = false;//not cograph
                    break;
                }
            }
        }
        assert(is_cograph);
    }
}  // namespace Koala
