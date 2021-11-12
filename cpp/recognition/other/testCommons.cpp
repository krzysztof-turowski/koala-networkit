#include "testCommons.h"

#include <chrono>
#include <cstdlib>
#include <iostream>

// #include "color.h"
#include "commons.h"
#include "oddHoles.h"
#include "perfect.h"

using namespace std::chrono;
using std::cerr;
using std::cout;
using std::default_random_engine;
using std::flush;
using std::make_pair;
using std::make_tuple;
using std::normal_distribution;

bool probTrue(double p) { return rand() / (RAND_MAX + 1.0) < p; }

void printGraph(const Graph &G) {
  cout << "   ";
  for (int i = 0; i < G.n; i++) {
    cout << i;
    if (i < 10) cout << " ";
  }
  cout << endl;
  for (int i = 0; i < G.n; i++) {
    cout << i;
    cout << ":";
    if (i < 10) cout << " ";
    for (int j = 0; j < G.n; j++) {
      cout << (G.areNeighbours(i, j) ? "X " : ". ");
    }
    cout << endl;
  }
  cout << endl;
}

Graph getRandomGraph(int size, double p) {
  vec<vec<int>> neighbours(size);
  for (int i = 0; i < size; i++) {
    for (int j = i + 1; j < size; j++) {
      if (probTrue(p)) {
        neighbours[i].push_back(j);
        neighbours[j].push_back(i);
      }
    }
  }

  return Graph(neighbours);
}

Graph getRandomPerfectGraph(int size, double p) {
  Graph G(0);
  do {
    G = getRandomGraph(size, p);
  } while (!isPerfectGraph(G));

  return G;
}

vec<Graph> getRandomGraphs(int size, double p, int howMany) {
  vec<Graph> ret;
  ret.reserve(howMany);
  for (int i = 0; i < howMany; i++) {
    ret.push_back(getRandomGraph(size, p));
  }

  return ret;
}

Graph getNonPerfectGraph(int holeSize, int reminderSize, double p) {
  if (holeSize < 5 || holeSize % 2 == 0) {
    throw invalid_argument("Trying to generate Non perfect graph with incorrect hole size.");
  }

  int size = holeSize + reminderSize;
  vec<vec<int>> neighbours(size);

  for (int i = 1; i < holeSize; i++) {
    neighbours[i - 1].push_back(i);
    neighbours[i].push_back(i - 1);
  }
  neighbours[0].push_back(holeSize - 1);
  neighbours[holeSize - 1].push_back(0);

  Graph reminder = getRandomGraph(reminderSize, p);
  for (int i = 0; i < reminderSize; i++) {
    for (int v : reminder[i]) {
      neighbours[i + holeSize].push_back(v + holeSize);
    }
  }

  for (int i = 0; i < holeSize; i++) {
    for (int j = holeSize; j < size; j++) {
      if (probTrue(p)) {
        neighbours[i].push_back(j);
        neighbours[j].push_back(i);
      }
    }
  }

  return Graph(neighbours).getShuffled();
}

Graph getBipariteGraph(int size, double p) {
  vec<vec<int>> neighbours(size);
  for (int i = 0; i < size / 2; i++) {
    for (int j = size / 2; j < size; j++) {
      if (probTrue(p)) {
        neighbours[i].push_back(j);
        neighbours[j].push_back(i);
      }
    }
  }

  return Graph(neighbours).getShuffled();
}

Graph getFullBinaryTree(int size) {
  vec<vec<int>> neighbors(size);
  for (int i = 1; i < size; i++) {
    neighbors[i].push_back(i / 2);
    neighbors[i / 2].push_back(i);
  }

  return Graph(neighbors).getShuffled();
}

Graph getGridWithMoves(int W, int H, const vec<int> &dx, const vec<int> &dy) {
  int n = W * H;

  assert(dx.size() == dy.size());

  vec<vec<int>> neighbors(n);
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
      for (int k = 0; k < dx.size(); k++) {
        int nx = x + dx[k];
        int ny = y + dy[k];
        if (nx >= 0 && nx < W && ny >= 0 && ny < H) {
          neighbors[x * H + y].push_back(nx * H + ny);
        }
      }
    }
  }

  return Graph(neighbors).getShuffled();
}

Graph getCityGrid(int W, int H) { return getGridWithMoves(W, H, {0, 0, 1, -1}, {1, -1, 0, 0}); }

Graph getKnightGraph(int W, int H) {
  vec<int> dx{1, 1, -1, -1, 2, 2, -2, -2};
  vec<int> dy{2, -2, 2, -2, 1, -1, 1, -1};

  return getGridWithMoves(W, H, dx, dy);
}

Graph getHypercube(int n) {
  vec<vec<int>> neighbors(n);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (__builtin_popcount(i ^ j) == 1) {
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
      }
    }
  }

  return Graph(neighbors).getShuffled();
}

Graph getRookGraph(int W, int H) {
  vec<int> dx;
  vec<int> dy;
  for (int i = 1; i <= W; i++) {
    dx.push_back(i);
    dy.push_back(0);

    dx.push_back(-i);
    dy.push_back(0);
  }

  for (int i = 1; i <= H; i++) {
    dy.push_back(i);
    dx.push_back(0);

    dy.push_back(-i);
    dx.push_back(0);
  }

  return getGridWithMoves(W, H, dx, dy);
}

Graph getSplitGraph(int n, double p) {
  vec<vec<int>> neighbors(n * 2);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) neighbors[i].push_back(j);
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = n; j < 2 * n; j++) {
      if (probTrue(p)) {
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
      }
    }
  }

  return Graph(neighbors).getShuffled();
}

void handler(int sig) {
  void *array[100];
  size_t size;

#ifdef __CYGWIN__
  fprintf(stderr, "Error: signal %d:\n", sig);
#else
  // get void*'s for all entries on the stack
  size = backtrace(array, 100);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  auto messages = backtrace_symbols(array, size);

  /* skip first stack frame (points here) */
  for (int i = 1; i < size && messages != NULL; ++i) {
    if ((messages[i][1] != 'l' || messages[i][2] != 'i' || messages[i][3] != 'b') &&
        (messages[i][1] != 'u' || messages[i][2] != 's' || messages[i][3] != 'r'))
      fprintf(stderr, "\t[bt]: (%d) %s\n", i, messages[i]);
  }
#endif

  exit(1);
}

void init(bool srandTime) {
  if (srandTime) {
    srand(time(NULL));
  }

  signal(SIGSEGV, handler);
  signal(SIGABRT, handler);
}

RaiiTimer::RaiiTimer(string msg) : msg(msg) {
  start_ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
}
RaiiTimer::~RaiiTimer() {
  nanoseconds end_ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
  double duration = (end_ns.count() - start_ns.count()) / 1e9;
  if (msg.size() > 0) cout << msg << ": " << duration << "s" << endl;
}
double RaiiTimer::getElapsedSeconds() {
  nanoseconds end_ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
  double duration = (end_ns.count() - start_ns.count()) / 1e9;
  return duration;
}

map<tuple<algos, bool, int>, double> sumTime;
// map<pair<int, bool>, double> sumClockTime;
map<tuple<algos, bool, int>, int> casesTested;
string algo_names[] = {"Perfect", "Naive", "CUDA Naive", "Cuda Perfect", "ALGO LAST", "CSDP Color"};

// map<pair<int, bool>, double> sumTimeNaive;
// map<pair<int, bool>, double> sumClockTimeNaive;
// map<pair<int, bool>, int> casesTestedNaive;

default_random_engine generator;
default_random_engine generatorWide;
normal_distribution<double> distribution(0.5, 0.05);
normal_distribution<double> distributionWide(0.5, 0.20);

/*bool testGraph(const Graph &G, vec<algos> algosToTest, vec<cuIsPerfectFunction> cuFunctions) {
  bool prevResult;
  bool result;

  for (int i = 0; i < algosToTest.size(); i++) {
    algos algo = algosToTest[i];
    StatsFactory::startTestCase(G, algo);

    switch (algo) {
      case algoPerfect:
        result = isPerfectGraph(G, true);
        break;

      case algoNaive:
        result = isPerfectGraphNaive(G, true);
        break;

      case algoCudaNaive:
        assert(cuFunctions[i] != nullptr);
        result = cuFunctions[i](G, true);
        break;

      case algoCudaPerfect:
        assert(cuFunctions[i] != nullptr);
        result = cuFunctions[i](G, true);
        break;

      default:
        throw invalid_argument("TestWithStats invalid argument");
    }

    if (i == 0) {
      prevResult = result;
    } else if (prevResult != result) {
      cerr << "Test Graph Error! " << algo_names[0] << " != " << algo_names[i] << endl;
      exit(1);
    }

    StatsFactory::endTestCase(result);
  }

  return result;
}*/

// void printStats() {
//   for (int algo = algoPerfect; algo < algo_last; algo++) {
//     int count = 0;
//     for (auto it = sumTime.begin(); it != sumTime.end(); it++) {
//       if (get<0>(it->first) == algo) count++;
//     }

//     if (count > 0) cout << algo_names[algo] << " recognition stats: " << endl;
//     for (auto it = sumTime.begin(); it != sumTime.end(); it++) {
//       if (get<0>(it->first) == algo) {
//         int cases = casesTested[it->first];
//         cout << "\tn=" << get<2>(it->first) << ", result=" << get<1>(it->first) << ", cases=" << cases
//              << ", avgTime=" << it->second / cases << endl;
//       }
//     }
//   }
// }

bool StatsFactory::curStarted = false;
int StatsFactory::testCaseNr = 0;
int StatsFactory::curTestPartNr = 0;
algos StatsFactory::curAlgo = algo_last;
int StatsFactory::curN = 0;
vec<string> StatsFactory::partNames = vec<string>();
vec<double> StatsFactory::curTime = vec<double>();
RaiiTimer StatsFactory::curTimer = RaiiTimer("");
RaiiTimer StatsFactory::curTimerOverall = RaiiTimer("");

map<string, int> StatsFactory::mapNameNr = map<string, int>();
map<int, string> StatsFactory::mapNrName = map<int, string>();
map<tuple<algos, bool, int, int>, int> StatsFactory::mapCount = map<tuple<algos, bool, int, int>, int>();
map<tuple<algos, bool, int, int>, double> StatsFactory::mapSumTime =
    map<tuple<algos, bool, int, int>, double>();

void StatsFactory::startTestCase(const Graph &G, algos algo_in) {
  testCaseNr++;
  curStarted = true;
  curAlgo = algo_in;
  curTestPartNr = 0;
  curN = G.n;
  curTimer = RaiiTimer("");
  curTimerOverall = RaiiTimer("");

  partNames.push_back("Overall");
  curTestPartNr++;
}

void StatsFactory::startTestCasePart(const string &name) {
  if (curTestPartNr > 0) {
    curTime.push_back(curTimer.getElapsedSeconds());
  }

  partNames.push_back(name);
  curTestPartNr++;

  curTimer = RaiiTimer("");
}

void StatsFactory::endTestCase(bool result) {
  double overallElapsed = curTimerOverall.getElapsedSeconds();

  if (curTestPartNr > 0) {
    curTime.push_back(curTimer.getElapsedSeconds());
  }

  for (int i = 0; i < partNames.size(); i++) {
    if (mapNameNr.count(partNames[i]) == 0) {
      mapNrName[mapNameNr.size()] = partNames[i];
      mapNameNr[partNames[i]] = mapNameNr.size();
    }
    auto t = make_tuple(curAlgo, result, curN, mapNameNr[partNames[i]]);
    // mapCount[t]++;
    mapSumTime[t] += curTime[i];
  }

  // This is so we sum multiple same parts in a algorithm
  set<string> setNames(partNames.begin(), partNames.end());
  for (string name : setNames) {
    auto t = make_tuple(curAlgo, result, curN, mapNameNr[name]);
    mapCount[t]++;
  }

  auto t = make_tuple(curAlgo, result, curN, mapNameNr["Overall"]);
  mapSumTime[t] += overallElapsed;

  curTime.clear();
  partNames.clear();
}

void StatsFactory::printStats2() {
  cout << "algorithm,\tn,\tresult,\tnum_runs,\t" << flush;
  for (int i = 0; i < mapNameNr.size(); i++) {
    cout << mapNrName[i] << flush;
    if (i + 1 < mapNrName.size()) cout << ",\t" << flush;
  }
  cout << endl;

  for (auto it = mapCount.begin(); it != mapCount.end(); it++) {
    auto t = it->first;
    int count = it->second;
    cout << algo_names[get<0>(t)] << ",\t" << get<2>(t) << ",\t" << get<1>(t) << ",\t" << count;

    for (int i = 0; i < mapNameNr.size(); i++) {
      auto t2 = t;
      get<3>(t2) = i;
      cout << ",\t";
      if (mapCount.count(t2) > 0) {
        cout << mapSumTime[t2] / mapCount[t2];
        if (mapCount[t2] != count) cout << "(" << mapCount[t2] << ")";
      } else {
        cout << " - ";
      }
    }

    for (; it != mapCount.end() && get<0>(it->first) == get<0>(t) && get<1>(it->first) == get<1>(t) &&
           get<2>(it->first) == get<2>(t);
         it++) {
      // cout << ",\t";
      // assert(count == mapCount[it->first]);
      // cout << mapSumTime[it->first] / count;
    }
    it--;

    cout << endl;
  }
}

double getDistr() {
  double distr = distribution(generator);
  if (distr <= 0 || distr > 1) distr = 0.5;

  return distr;
}

double getDistrWide() {
  double distr = distributionWide(generatorWide);
  if (distr <= 0 || distr > 1) distr = 0.5;

  return distr;
}

void printTimeHumanReadable(double time, bool use_cerr) {
  auto &out = use_cerr ? cerr : cout;

  double s = time;
  int h = s / (60 * 60);
  s -= h * (60 * 60);
  if (h != 0) {
    out << h << "h";
  }

  int m = s / 60;
  s -= m * 60;
  if (m != 0) {
    out << m << "m";
  }

  out << static_cast<int>(s) + 1 << "s";
}

RaiiProgressBar::RaiiProgressBar(int allTests, bool use_cerr) : allTests(allTests), use_cerr(use_cerr) {
  start_ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
  update(0);
}

RaiiProgressBar::~RaiiProgressBar() {
  if (use_cerr)
    cerr << endl;
  else
    cout << endl;
}

int RaiiProgressBar::getFilled(int testsDone) {
  double progress = testsDone / (static_cast<double>(allTests));
  return width * progress;
}

void RaiiProgressBar::update(int testsDone) {
  auto &out = use_cerr ? cerr : cout;

  int toFill = getFilled(testsDone);
  if (testsDone < 10 || testsDone == allTests || toFill != getFilled(testsDone - 1)) {
    out << "[";
    for (int i = 0; i < width; i++) {
      out << (i < toFill ? "X" : " ");
    }
    out << "]";
    if (testsDone > 0) {
      out << " (about ";
      nanoseconds end_ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
      double timeElapsed = (end_ns.count() - start_ns.count()) / 1e9;
      double timeRemaining = timeElapsed * (allTests - testsDone) / testsDone;
      printTimeHumanReadable(timeRemaining, use_cerr);
      out << " left)";
    }
    out << "\r" << flush;
  }
}
