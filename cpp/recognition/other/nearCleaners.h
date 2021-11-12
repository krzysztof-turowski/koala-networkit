#pragma once

#include <boost/dynamic_bitset.hpp>
#include <set>
#include "commons.h"

// Returns true if G has an odd hole.
// Returns false if there is no shortest odd hole C such that X is a near-cleaner for C.
// G should contain no pyramid or jewel.
bool containsOddHoleWithNearCleanerX(const Graph &G, const set<int> &sX, const vec<vec<int>> &triplePaths,
                                     bool gatherStats = false);

bool isRelevantTriple(const Graph &G, int a, int b, int c);

// Returns X(a, b, c).
// (a, b, c) should be relevant triple.
boost::dynamic_bitset<ul> getXforRelevantTriple(const Graph &G, const Graph &GC, int a, int b, int c);

// 9.2
set<boost::dynamic_bitset<ul>> getPossibleNearCleaners(const Graph &G, const Graph &GC,
                                                       bool gatherStats = false);
