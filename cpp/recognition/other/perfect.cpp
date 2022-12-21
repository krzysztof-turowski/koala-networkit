#include "perfect.h"
#include "commons.h"
#include "nearCleaners.h"
#include "oddHoles.h"
#include "pyramids.h"
#include "testCommons.h"
#include <networkit/graph/Graph.hpp>

Koala::PerfectGraphRecognition::State containsSimpleProhibited(const Graph &G, bool gatherStats) {
  if (gatherStats) StatsFactory::startTestCasePart("T2");
  if (containsT2(G)) return Koala::PerfectGraphRecognition::State::HAS_T2;

  if (gatherStats) StatsFactory::startTestCasePart("T3");
  if (containsT3(G)) return Koala::PerfectGraphRecognition::State::HAS_T3;

  return Koala::PerfectGraphRecognition::State::UNKNOWN;
}

Koala::PerfectGraphRecognition::State checkPerfectGraph(const NetworKit::Graph &graph, bool gatherStats) {
  Graph G(graph);
  return checkPerfectGraph(G, gatherStats);
}

Koala::PerfectGraphRecognition::State checkPerfectGraph(const Graph &G, bool gatherStats) {
  const bool printInterestingGraphs = false;

  if (gatherStats) StatsFactory::startTestCasePart("Simple Structures");

  auto simple_test = containsSimpleProhibited(G, gatherStats);
  if (simple_test != Koala::PerfectGraphRecognition::State::UNKNOWN) {
    return simple_test;
  }

  Graph GC = G.getComplement();
  auto simple_test_complement = containsSimpleProhibited(GC, gatherStats);
  if (simple_test_complement != Koala::PerfectGraphRecognition::State::UNKNOWN) {
    return simple_test_complement;
  }

  if (gatherStats) StatsFactory::startTestCasePart("Get Near Cleaners");
  auto Xs = getPossibleNearCleaners(G, GC, false);

  if (gatherStats) StatsFactory::startTestCasePart("Test NC Rest");
  vec<vec<int>> triplePaths = getAllPaths(G, 3);

  for (auto X : Xs) {
    if (containsOddHoleWithNearCleanerX(G, bitsetToSet(X), triplePaths, gatherStats)) {
      if (printInterestingGraphs) cout << "Interesting graph: " << G << endl;
      return Koala::PerfectGraphRecognition::State::HAS_NEAR_CLEANER_ODD_HOLE;
    }
  }

  if (gatherStats) StatsFactory::startTestCasePart("Get Near Cleaners");
  auto XsC = getPossibleNearCleaners(GC, G, false);

  if (gatherStats) StatsFactory::startTestCasePart("Test NC Rest");
  vec<vec<int>> CTriplePaths = getAllPaths(GC, 3);

  for (auto X : XsC) {
    if (containsOddHoleWithNearCleanerX(GC, bitsetToSet(X), CTriplePaths, gatherStats)) {
      if (printInterestingGraphs) cout << "Interesting graph: " << G << endl;
      return Koala::PerfectGraphRecognition::State::HAS_NEAR_CLEANER_ODD_HOLE;
    }
  }

  return Koala::PerfectGraphRecognition::State::PERFECT;
}

bool isPerfectGraph(const NetworKit::Graph &graph, bool gatherStats) {
  return checkPerfectGraph(graph, gatherStats) == Koala::PerfectGraphRecognition::State::PERFECT;
}

bool isPerfectGraph(const Graph &G, bool gatherStats) {
  return checkPerfectGraph(G, gatherStats) == Koala::PerfectGraphRecognition::State::PERFECT;
}

bool isPerfectGraphNaive(const NetworKit::Graph &graph, bool gatherStats) {
  Graph G(graph);
  return isPerfectGraphNaive(G, gatherStats);
}

bool isPerfectGraphNaive(const Graph &G, bool gatherStats) {
  return !containsOddHoleNaive(G, gatherStats) && !containsOddHoleNaive(G.getComplement(), false);
}
