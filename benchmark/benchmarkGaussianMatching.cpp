#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>

#include <io/DimacsGraphReader.hpp>

#include <matching/MaximumMatching.hpp>
#include <matching/gaussian_matching/BipartiteGaussianMatching.hpp>
#include <matching/gaussian_matching/GeneralGaussianMatching.hpp>
#include <matching/gaussian_matching/NaiveGaussianMatching.hpp>

using namespace std;
using namespace chrono;
using namespace NetworKit;

struct MatchingBenchmark {
  MatchingBenchmark(string name, int matchingSize, int perfectSize,
                    steady_clock::time_point begin,
                    steady_clock::time_point end)
      : name(name), matchingSize(matchingSize), perfectSize(perfectSize),
        duration(duration_cast<microseconds>(end - begin).count()) {}

  string name;
  int matchingSize, perfectSize;
  uint64_t duration;

  friend ostream &operator<<(ostream &os, const MatchingBenchmark &mb) {
    os << "Benchmark for " << mb.name << "algorithm:\n";
    os << "\tMatching size:\t" << mb.matchingSize << "/" << mb.perfectSize
       << '\n';
    os << "\tTime:\t" << mb.duration << "μs\n";
    return os;
  }
};

MatchingBenchmark runGeneralGaussianMatching(Graph &G) {
  auto begin = steady_clock::now();
  Koala::GeneralGaussianMatching algo(G);
  algo.run();
  auto matchingSize = algo.getMatching().size();  // only one of (u,v), (v,u)
  auto end = steady_clock::now();
  return {"GeneralGaussianMatching", matchingSize, G.numberOfNodes() / 2, begin,
          end};
}

MatchingBenchmark runBipartiteGaussianMatching(Graph &G) {
  auto begin = steady_clock::now();
  Koala::GeneralGaussianMatching algo(G);
  algo.run();
  auto matchingSize = algo.getMatching().size();  // only one of (u,v), (v,u)
  auto end = steady_clock::now();
  return {"BipartiteGaussianMatching", matchingSize, G.numberOfNodes() / 2,
          begin, end};
}

MatchingBenchmark runNaiveGaussianMatching(Graph &G) {
  auto begin = steady_clock::now();
  Koala::NaiveGaussianMatching algo(G);
  algo.run();
  auto matchingSize = algo.getMatching().size() / 2;  // both (u,v), (v,u)
  auto end = steady_clock::now();
  return {"NaiveGaussianMatching", matchingSize, G.numberOfNodes() / 2, begin,
          end};
}

MatchingBenchmark runEndmondsMatching(Graph &G) {
  auto begin = steady_clock::now();
  Koala::EdmondsMaximumMatching algo(G, true,
                                     Koala::BlossomMaximumMatching::empty);
  algo.run();
  auto matchingSize = algo.getMatching().size() / 2;  // both (u,v), (v,u)
  auto end = steady_clock::now();
  return {"EdmondsMaximumMatching", matchingSize, G.numberOfNodes() / 2, begin,
          end};
}

MatchingBenchmark runGabowMatching(Graph &G) {
  auto begin = steady_clock::now();
  Koala::GabowMaximumMatching algo(G, true,
                                   Koala::BlossomMaximumMatching::empty);
  algo.run();
  auto matchingSize = algo.getMatching().size() / 2;  // both (u,v), (v,u)
  auto end = steady_clock::now();
  return {"GabowMaximumMatching", matchingSize, G.numberOfNodes() / 2, begin,
          end};
}

Graph readGraph(const string &filepath) {
  return Koala::DimacsGraphReader().read(filepath);
}

void benchmarkAverage(function<MatchingBenchmark(Graph &)> runAlgorithm,
                      vector<vector<string>> testCases) {
  vector<uint64_t> benchmark;
  string name;
  for (int i = 0; i < testCases.size(); ++i) {
    auto testCase = testCases[i];
    uint64_t totalTime = 0;
    cerr << "Running tc " << i << "...\n";
    for (auto testPath : testCase) {
      auto G = readGraph(testPath);
      G.indexEdges();
      cerr << "\tRunning " << testPath << "...";
      auto result = runAlgorithm(G);
      cerr << "\tdone\n";
      totalTime += result.duration;
      name = result.name;
    }
    cerr << "done\n";
    benchmark.push_back(totalTime / testCase.size());
  }

  cout << name << '\n';
  for (auto time : benchmark) {
    cout << '\t' << time << "μs\n";
  }
}

int main() {
  string types[] = {"gen", "bp"};
  pair<int, int> sizes[] = {{100, 2000}, {1000, 200000}};
  int K = 10;

  vector<vector<string>> testCases;
  vector<vector<string>> bpTestCases;
  for (auto type : types) {
    for (auto [N, M] : sizes) {
      vector<string> testCase;
      for (int i = 0; i < K; ++i) {
        std::stringstream ss;
        ss << "./input/test_" << N << "_" << M << "_" << i << "." << type;
        testCase.push_back(ss.str());
      }
      testCases.push_back(testCase);
      if (type == "bp") {
        bpTestCases.push_back(testCase);
      }
    }
  }

  benchmarkAverage(runEndmondsMatching, testCases);
  benchmarkAverage(runGabowMatching, testCases);
  benchmarkAverage(runGeneralGaussianMatching, testCases);
  benchmarkAverage(runBipartiteGaussianMatching, bpTestCases);
  benchmarkAverage(runNaiveGaussianMatching, testCases);

  return 0;
}
