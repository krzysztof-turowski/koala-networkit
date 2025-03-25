#pragma once

#include <vector>

namespace Koala {

class FenwickTree {
 public:
    int sum(int index) {
        index++;
        int res = 0;
        for (; index > 0; index -= index & (-index))
            res += T[index];
        return res;
    }

    void add(int index, int val) {
        index++;
        for (; index < T.size(); index += index & (-index))
            T[index] += val;
    }

    void reset(int size) {
        T.resize(size + 1);
        for (int i = 0; i <= size; ++i)
            T[i] = 0;
    }

 private:
    std::vector<int> T;
};

} /* namespace Koala */
