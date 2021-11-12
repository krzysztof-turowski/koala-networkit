#pragma once

#include <tuple>
#include "commons.h"

bool checkPrerequisites(const Graph &G, const vec<int> &b, const int a, const vec<int> &s);

bool vectorsCutEmpty(vec<int>::iterator aBegin, vec<int>::iterator aEnd, vec<int>::iterator bBegin,
                     vec<int>::iterator bEnd);

bool noEdgesBetweenVectors(const Graph &G, vec<int>::iterator aBegin, vec<int>::iterator aEnd,
                           vec<int>::iterator bBegin, vec<int>::iterator bEnd);

// If pyramide is found, returns [b1, b2, b3], a, [P1, P2, P3]
// Returns empty vector, -1 and emty vector if none found
tuple<vec<int>, int, vec<vec<int>>> findPyramid(const Graph &G);

// Checks if [b1, b2, b3], a, [P1, P2, P3] is a pyramide
// used only in tests
bool isPyramid(const Graph &G, vec<int> b, int a, vec<vec<int>> P);

bool containsPyramid(const Graph &G);
