

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
        if(x){
            if(x == 1)return CographRecognition::State::COND_1;
            else if(x == 2)return CographRecognition::State::COND_2;
            else if(x == 3)return CographRecognition::State::COND_3;
            else if(x == 4)return CographRecognition::State::COND_4;
            else if(x == 5)return CographRecognition::State::COND_5;
            else if(x == 6)return CographRecognition::State::COND_6;
        }
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




}  // namespace Koala
