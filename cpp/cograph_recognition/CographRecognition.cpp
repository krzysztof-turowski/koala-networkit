#include <list>
#include <graph/GraphTools.hpp>
#include <cograph_recognition/CographRecognition.hpp>

namespace Koala {

CographRecognition::CographRecognition(NetworKit::Graph &graph)
    : graph(graph), is_complement_reducible(State::UNKNOWN) {

     }
bool CographRecognition::isComplementReducible() const {
    assureFinished();
    return is_complement_reducible == State::COMPLEMENT_REDUCIBLE;
}

CographRecognition::State CographRecognition::getState() const {
    assureFinished();
    return is_complement_reducible;
}
    CographRecognition::State recognition(NetworKit::Graph &graph){
        return  CographRecognition::Cograph_Recognition(graph);
    }
void CographRecognition::run() {
    hasRun = true;
    if (is_complement_reducible != State::UNKNOWN) {
        return;
    }
    is_complement_reducible = recognition(graph);
    if (is_complement_reducible != State::UNKNOWN) {
            return;
    }

    is_complement_reducible = State::COMPLEMENT_REDUCIBLE;
}

    void CographRecognition::check() const {
        assureFinished();
        bool iscograph = true;
        for(auto e1 : graph.edgeRange()){
            for(auto e2 : graph.edgeRange()){
                auto x = e1.u;
                auto y = e1.v;
                auto u = e2.u;
                auto v = e2.v;
                if(!(x != u && x != v && y != u && y != v))continue;
                if(graph.hasEdge(y, u) && !graph.hasEdge(x, u) && !graph.hasEdge(x, v) && !graph.hasEdge(y, v)) {
                    iscograph = false;//not cograph
                    break;
                }
            }
        }
        assert(iscograph);
    }


}  // namespace Koala
