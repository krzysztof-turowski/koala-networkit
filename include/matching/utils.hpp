#pragma once

#include <cmath>
#include <random>

namespace Koala {
    constexpr double EPS = 1e-8;

    inline bool eq(double a, double b) {
        return fabs(a - b) <= EPS;
    }

    inline double generateRandom() {
        return (double)rand() / RAND_MAX; //TODO
    }
}