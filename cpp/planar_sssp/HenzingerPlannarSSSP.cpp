#include "planar_sssp/HenzingerPlanarSSSP.hpp"
#include "planar_sssp/PlanarUtils.hpp"

#include <set>
#include <unordered_map>
#include <utility>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using std::cout;
using std::endl;

using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
    void HenzingerPlanarSSSP::run() {

        hasRun = true;
        return;
    }
}