//
// Created by milana on 19.03.24.
//
#include <cassert>
#include <iostream>
#include <map>

#include "io/G6GraphReader.hpp"
#include "cograph_recognition/CographRecognition.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;
int main() {
    std::ifstream fin("../../input/cographConnected14.g6");
    std::map<Koala::CographRecognition::State, int> classification;
    std::string types[] = {
            "UNKNOWN",
            "COMPLEMENT_REDUCIBLE",
            "COND_1",
            "COND_2",
            "COND_3",
            "COND_4",
            "COND_5",
            "COND_6"
    };
    cout<<"why?"<<endl;
    std::cout << "Current path is " << fs::current_path() <<" "<<fin.is_open()<<std::endl;
    int cnt = 0;
    while (cnt !=2 ) {
        cnt++;
        std::string line;
        fin >> line;
        if (!std::cin.good()) {
            break;
        }
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::CographRecognition(G);
        cout<<"aaa"<<endl;
        recognize.run();
        cout<<"bbb"<<endl;
        classification[recognize.getState()]++;
        cout<<"ccc"<<endl;
        recognize.check();
        cout<<"ddd"<<endl;
    }
    for (const auto &[k, v] : classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return 0;
}
