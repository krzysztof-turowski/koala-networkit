#pragma once

#include "commons.h"

#include <recognition/PerfectGraphRecognition.hpp>

Koala::PerfectGraphRecognition::State checkPerfectGraph(const NetworKit::Graph &graph, bool gatherStats = false);
Koala::PerfectGraphRecognition::State checkPerfectGraph(const Graph &G, bool gatherStats = false);

bool isPerfectGraph(const NetworKit::Graph &graph, bool gatherStats = false);
bool isPerfectGraph(const Graph &G, bool gatherStats = false);

Koala::PerfectGraphRecognition::State containsSimpleProhibited(const Graph &G, bool gatherStats = false);

bool isPerfectGraphNaive(const NetworKit::Graph &graph, bool gatherStats = false);
bool isPerfectGraphNaive(const Graph &G, bool gatherStats = false);
