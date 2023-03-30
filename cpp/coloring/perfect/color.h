#include <tuple>
#include "commons.h"

// Returns: n, m, from, to suitable for theta function from csdp. Nodes start from 1, unlike in Graph.
// isNodeRemoved is a vector of booleans. If isNodeRemoved[i] is present and true, this node is removed from
// graph before processing
tuple<int, int, vec<int>, vec<int>> getGraphEdges(const Graph &G, const vec<int> &isNodeRemoved = vec<int>());

int getTheta(const Graph &G, const vec<int> &isNodeRemoved = vec<int>(), bool gatherStats = false);
int getOmega(const Graph &G, bool gatherStats = false);

bool isStableSet(const Graph &G, vec<int> nodes);

// Returned vector is sorted
vec<int> getMaxCardStableSet(const Graph &G, bool gatherStats = false);
// Returned vector is sorted
vec<int> getMaxCardClique(const Graph &G, bool gatherStats = false);

// Returns a Stable Set intersecting all given max-card Cliques K[0], ... K.back()
// Returned vector is sorted.
vec<int> getSSIntersectingCliques(const Graph &G, vec<vec<int>> K, bool gatherStats = false);

// Returns a Stable Set intersecting all max card Cliques
vec<int> getSSIntersectingAllMaxCardCliques(const Graph &G, bool gatherStats = false);

vec<int> color(const Graph &G, bool gatherStats = false);
