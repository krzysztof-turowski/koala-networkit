#pragma once

#include <networkit/Globals.hpp>
#include <vector>
namespace Koala {
class LiptonTarjanEdgeScanner {
  public:
    LiptonTarjanEdgeScanner(NetworKit::node beg, NetworKit::node prev, int rotDir,
                            std::vector<NetworKit::node> &cyclePrev);

  private:
    NetworKit::node beg, cur = beg, prev;
    int embCur = -1, embStop = -1;
    int rotDir;
    double cost = 0.0;
    bool forward = false, done = false, isInitRequired = true;
};
} // namespace Koala
