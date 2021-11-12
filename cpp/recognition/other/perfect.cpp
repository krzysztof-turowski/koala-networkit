#include "perfect.h"
#include "commons.h"
#include "jewels.h"
#include "nearCleaners.h"
#include "oddHoles.h"
#include "pyramids.h"
#include "testCommons.h"
#include <networkit/graph/Graph.hpp>

bool containsSimpleProhibited(const Graph &G, bool gatherStats) {
  if (gatherStats) StatsFactory::startTestCasePart("Jewel");
  if (containsJewelNaive(G)) return true;

  if (gatherStats) StatsFactory::startTestCasePart("Pyramid");
  if (containsPyramid(G)) return true;

  if (gatherStats) StatsFactory::startTestCasePart("T1");
  if (containsT1(G)) return true;

  if (gatherStats) StatsFactory::startTestCasePart("T2");
  if (containsT2(G)) return true;

  if (gatherStats) StatsFactory::startTestCasePart("T3");
  if (containsT3(G)) return true;

  return false;
}

bool isPerfectGraph(const NetworKit::Graph &graph, bool gatherStats) {
  Graph G(graph);
  return isPerfectGraph(G, gatherStats);
}

bool isPerfectGraph(const Graph &G, bool gatherStats) {
  const bool printInterestingGraphs = false;

  if (gatherStats) StatsFactory::startTestCasePart("Simple Structures");

  Graph GC = G.getComplement();
  if (containsSimpleProhibited(G, gatherStats) || containsSimpleProhibited(GC, gatherStats)) return false;

  if (gatherStats) StatsFactory::startTestCasePart("Get Near Cleaners");
  auto Xs = getPossibleNearCleaners(G, GC, false);

  if (gatherStats) StatsFactory::startTestCasePart("Test NC Rest");
  vec<vec<int>> triplePaths = getAllPaths(G, 3);

  for (auto X : Xs) {
    if (containsOddHoleWithNearCleanerX(G, bitsetToSet(X), triplePaths, gatherStats)) {
      if (printInterestingGraphs) cout << "Interesting graph: " << G << endl;
      return false;
    }
  }

  if (gatherStats) StatsFactory::startTestCasePart("Get Near Cleaners");
  auto XsC = getPossibleNearCleaners(GC, G, false);

  if (gatherStats) StatsFactory::startTestCasePart("Test NC Rest");
  vec<vec<int>> CTriplePaths = getAllPaths(GC, 3);

  for (auto X : XsC) {
    if (containsOddHoleWithNearCleanerX(GC, bitsetToSet(X), CTriplePaths, gatherStats)) {
      if (printInterestingGraphs) cout << "Interesting graph: " << G << endl;
      return false;
    }
  }

  return true;
}

bool isPerfectGraphNaive(const NetworKit::Graph &graph, bool gatherStats) {
  Graph G(graph);
  return isPerfectGraphNaive(G, gatherStats);
}

bool isPerfectGraphNaive(const Graph &G, bool gatherStats) {
  return !containsOddHoleNaive(G, gatherStats) && !containsOddHoleNaive(G.getComplement(), false);
}
