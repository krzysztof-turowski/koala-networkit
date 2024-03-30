//
// Created by milana on 19.03.24.
//
#include <iostream>
#include <map>
#include "io/G6GraphReader.hpp"
#include "cograph_recognition/CographRecognition.hpp"
#include <fstream>
using namespace std;
std::string types[] = {
        "UNKNOWN",
        "COGRAPH",
        "COND_1",
        "COND_2",
        "COND_3",
        "COND_4",
        "COND_5",
        "COND_6"
};
pair<int,int> test(int j, string s, string s1){
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
    }
    cout <<"size: "<<j<< " total: "<<cnt<< endl;
    for (const auto &[k, v]: classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
    return {cnt,classification[Koala::CographRecognition::State::COGRAPH]};
}
void auto_test(){
    for(int j = 3; j <= 9; j++){

        pair<int,int> cnt = test(j, "../../input/cographConnected", ".g6");
        pair<int,int> cnt1 = test(j, "../../input/graph","c.g6");
        if(cnt.first == cnt.second && cnt1.second == cnt.first)cout<<"OK"<<endl;
        else cout<<"FAIL"<<endl;
    }
    for(int j = 3; j <= 15; j++) {
        test(j,"../../input/cographConnected", ".g6");
    }
}
void manual_test(){
    std::map<Koala::CographRecognition::State, int> classification;
    int cnt = 0;
    while (true) {
        std::string line;
        cin >> line;
        if (!cin.good()) {
            break;
        }
        cnt++;
        NetworKit::Graph G = Koala::G6GraphReader().readline(line);
        auto recognize = Koala::CographRecognition(G);
        recognize.run();
        classification[recognize.getState()]++;
    }
    cout << " total: "<<cnt<< endl;
    for (const auto &[k, v]: classification) {
        std::cout << types[static_cast<int>(k)] << ": " << v << std::endl;
    }
}
int main() {
    ///auto_test();//to check the correctness of algorithm
    manual_test();
    return 0;
}
