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
    for(int j = 3; j <= 10; j++){
        string s = "../../input/cographConnected", s1 = ".g6";
        string s2 = "";
        if(j >= 10)s2 += "1";
        s2 += '0' + j % 10;
        std::ifstream fin(s + s2 + s1);
        std::map<Koala::CographRecognition::State, int> classification;
        int cnt = 0;
        while (true) {


            std::string line;
            fin >> line;
            if (!fin.good()) {
                break;
            }
            cnt++;
            NetworKit::Graph G = Koala::G6GraphReader().readline(line);
            auto recognize = Koala::CographRecognition(G);
            recognize.run();
            classification[recognize.getState()]++;
            recognize.check();
            if (recognize.getState() != Koala::CographRecognition::State::COMPLEMENT_REDUCIBLE)cout << cnt << endl;
        }
        cout <<"size: "<<j<< " total: "<<cnt<< endl;
        for (const auto &[k, v]: classification) {
            std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
        }



        s = "../../input/graph", s1 = "c.g6";
        s2 = "";
        if(j >= 10)s2 += "1";
        s2 += '0' + j % 10;
        std::ifstream din(s + s2 + s1);
        std::map<Koala::CographRecognition::State, int> classification1;
        int cnt1 = 0;
        while (true) {


            std::string line;
            din >> line;
            if (!din.good()) {
                break;
            }
            cnt1++;
            NetworKit::Graph G = Koala::G6GraphReader().readline(line);
            auto recognize = Koala::CographRecognition(G);
            recognize.run();
            classification1[recognize.getState()]++;
            recognize.check();
        }
        cout <<"size: "<<j<< " total: "<<cnt1<< endl;
        for (const auto &[k, v]: classification1) {
            std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
        }
        if(cnt == classification[Koala::CographRecognition::State::COMPLEMENT_REDUCIBLE] && classification1[Koala::CographRecognition::State::COMPLEMENT_REDUCIBLE] == cnt)cout<<"OK"<<endl;
        else cout<<"FAIL"<<endl;
    }
    for(int j = 3; j <= 14; j++) {
        string s = "../../input/cographConnected", s1 = ".g6";
        string s2 = "";
        if(j >= 10)s2 += "1";
        s2 += '0' + j % 10;
        std::ifstream fin(s + s2 + s1);
        std::map<Koala::CographRecognition::State, int> classification;
        int cnt = 0;
        while (true) {


            std::string line;
            fin >> line;
            if (!fin.good()) {
                break;
            }
            cnt++;
            NetworKit::Graph G = Koala::G6GraphReader().readline(line);
            auto recognize = Koala::CographRecognition(G);
            recognize.run();
            classification[recognize.getState()]++;
            recognize.check();
            if (recognize.getState() != Koala::CographRecognition::State::COMPLEMENT_REDUCIBLE)cout << cnt << endl;
        }
        cout <<"size: "<<j<< " total: "<<cnt<< endl;
        for (const auto &[k, v]: classification) {
            std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
        }
    }
    return 0;
}
