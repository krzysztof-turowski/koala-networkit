

#include <list>
#include <set>

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
        int x = CographRecognition::Cograph_Recognition(graph);
        std::cout<<x<<std::endl;
        if(x){
            if(x == 1)return CographRecognition::State::COND_1;
            else if(x == 2)return CographRecognition::State::COND_2;
            else if(x == 3)return CographRecognition::State::COND_3;
            else if(x == 4)return CographRecognition::State::COND_4;
            else if(x == 5)return CographRecognition::State::COND_5;
            else if(x == 6)return CographRecognition::State::COND_6;
        }
        return CographRecognition::State::UNKNOWN;
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
                if(graph.hasEdge(y, u) && !graph.hasEdge(x, u) && !graph.hasEdge(x, v) && !graph.hasEdge(y, v)) {
                    iscograph = false;//not cograph
                    break;
                }
            }
        }
        assert(iscograph);
    }

}  // namespace Koala
