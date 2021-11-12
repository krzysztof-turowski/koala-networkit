#pragma once

#include "commons.h"

bool isPerfectGraph(const NetworKit::Graph &graph, bool gatherStats = false);
bool isPerfectGraph(const Graph &G, bool gatherStats = false);

bool containsSimpleProhibited(const Graph &G, bool gatherStats = false);

bool isPerfectGraphNaive(const NetworKit::Graph &graph, bool gatherStats = false);
bool isPerfectGraphNaive(const Graph &G, bool gatherStats = false);
