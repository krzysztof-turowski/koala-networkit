#pragma once

#include "commons.h"

bool isJewel(const Graph &G, const vec<int> &v);

// returns [v1, ..., v5] or empty vector if none found
// returned[0] = v1, returned[1] = v2, ...
// Runv in N^7, instead of N^6 in paper, but much simpler
vec<int> findJewelNaive(const Graph &G);
bool containsJewelNaive(const Graph &G);

// bool containsJewel(const Graph &G);
