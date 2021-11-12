#pragma once

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <map>
#include <random>
#include <string>
#include "commons.h"

#ifndef __CYGWIN__
#include <execinfo.h>
#endif

using std::get;
using std::invalid_argument;
using std::map;
using namespace std::chrono;

// Whether to run big tests. These take more time.
// const bool bigTests = false;  // TODO(Adrian) make bigTests big again (perf)
const bool bigTests = true;
const int _max_threads_to_run = 10000;

bool probTrue(double p);
void printGraph(const Graph &G);
Graph getRandomGraph(int size, double p);
// will loop when impossible or unlikely.
Graph getRandomPerfectGraph(int size, double p);
vec<Graph> getRandomGraphs(int size, double p, int howMany);

Graph getNonPerfectGraph(int holeSize, int reminderSize, double p);

// Returns biparite graph, two equal layers (+-1), each edge between layers has probability of p.
Graph getBipariteGraph(int size, double p);

Graph getFullBinaryTree(int size);
Graph getGridWithMoves(int W, int H, const vec<int> &dx, const vec<int> &dy);
Graph getCityGrid(int W, int H);
Graph getKnightGraph(int W, int H);
Graph getHypercube(int k);
Graph getRookGraph(int W, int H);
Graph getSplitGraph(int n, double p);

void handler(int sig);
void init(bool srandTime = false);

struct RaiiTimer {
  explicit RaiiTimer(string msg = "");
  ~RaiiTimer();

  double getElapsedSeconds();

 private:
  nanoseconds start_ns;
  string msg;
};

enum algos { algoPerfect, algoNaive, algoCudaNaive, algoCudaPerfect, algo_last, algo_color};
extern string algo_names[];

typedef bool (*cuIsPerfectFunction)(const Graph &G, bool gatherStats);
//bool testGraph(const Graph &G, vec<algos> algosToTest, vec<cuIsPerfectFunction> cuFunctions = {});

struct StatsFactory {
  static void startTestCase(const Graph &G, algos algo);
  static void startTestCasePart(const string &name);
  static void endTestCase(bool result);
  static void printStats2();

 private:
  static bool curStarted;
  static int testCaseNr;
  static int curTestPartNr;
  static algos curAlgo;
  static int curN;
  static vec<string> partNames;
  static vec<double> curTime;
  static RaiiTimer curTimer;
  static RaiiTimer curTimerOverall;

  static map<string, int> mapNameNr;
  static map<int, string> mapNrName;
  static map<std::tuple<algos, bool, int, int>, int> mapCount;
  static map<std::tuple<algos, bool, int, int>, double> mapSumTime;
};

double getDistr();
double getDistrWide();

void printTimeHumanReadable(int64_t time, bool use_cerr = true);

struct RaiiProgressBar {
  explicit RaiiProgressBar(int allTests, bool use_cerr = true);
  ~RaiiProgressBar();

  void update(int testsDone);

 private:
  int allTests;
  bool use_cerr;
  nanoseconds start_ns;
  const int width = 180;

  int getFilled(int testsDone);
};
