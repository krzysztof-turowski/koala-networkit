/*
 * RankTwoRelaxationMaxCut.cpp
 *
 * Solution for the Max-Cut problem using Rank Two Relaxation method.
 * Created on: 26.03.2024
 * Author: Michał Miziołek
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>

#include <max_cut/RankTwoRelaxationMaxCut.hpp>

namespace Koala {

RankTwoRelaxationMaxCut::RankTwoRelaxationMaxCut(const std::vector<std::vector<int>>& graphInput)
    : graph(graphInput), numberOfVertices(graphInput.size()), bestCutValue(std::numeric_limits<int>::min()) {
    theta.resize(numberOfVertices);
    distributeThetaEvenly();
}

void RankTwoRelaxationMaxCut::distributeThetaEvenly() {
    for (int i = 0; i < numberOfVertices; ++i) {
        theta[i] = 2 * M_PI * i / numberOfVertices;
    }
}

double RankTwoRelaxationMaxCut::computeCutValue(const std::vector<int>& x) {
    double cutValue = 0;
    for (int i = 0; i < numberOfVertices; ++i) {
        for (int j = 0; j < numberOfVertices; ++j) {
            if (x[i] != x[j]) {
                cutValue += graph[i][j];
            }
        }
    }
    return cutValue / 2;
}

std::vector<int> RankTwoRelaxationMaxCut::procedureCut() {
    double bestValue = -std::numeric_limits<double>::infinity();
    std::vector<int> bestCut(numberOfVertices), x(numberOfVertices);

    for (double alpha = 0; alpha <= M_PI; alpha += 0.01) {
        for (int i = 0; i < numberOfVertices; ++i) {
            x[i] = (theta[i] >= alpha && theta[i] < alpha + M_PI) ? 1 : -1;
        }

        double value = computeCutValue(x);
        if (value > bestValue) {
            bestValue = value;
            bestCut = x;
        }
    }
    return bestCut;
}

std::vector<double> RankTwoRelaxationMaxCut::calculateGradient(const std::vector<double>& theta) {
    std::vector<double> gradient(numberOfVertices, 0.0);
    for (int j = 0; j < numberOfVertices; ++j) {
        for (int k = 0; k < numberOfVertices; ++k) {
            gradient[j] += graph[k][j] * sin(theta[k] - theta[j]);
        }
    }
    return gradient;
}

void RankTwoRelaxationMaxCut::gradientDescent(double alpha, int maxIterations) {
    for (int iter = 0; iter < maxIterations; ++iter) {
        std::vector<double> gradient = calculateGradient(theta);
        for (int j = 0; j < numberOfVertices; ++j) {
            theta[j] -= alpha * gradient[j];
        }
    }

    for (int i = 0; i < numberOfVertices; ++i) {
        if (theta[i] < 0) theta[i] += 2 * M_PI;
        if (theta[i] >= 2 * M_PI) theta[i] -= 2 * M_PI;
    }

    std::sort(theta.begin(), theta.end());
}

void RankTwoRelaxationMaxCut::perturbTheta() {
    distributeThetaEvenly();
    gradientDescent(0.001, 1000000);
}

void RankTwoRelaxationMaxCut::solve() {
    perturbTheta();
    bestSet = procedureCut();
    bestCutValue = computeCutValue(bestSet);

    if (bestSet[numberOfVertices - 1] != -1) {
        bestSet[numberOfVertices - 1] = -1;
        double candidateValue = computeCutValue(bestSet);
        if (candidateValue > bestCutValue) {
            bestCutValue = candidateValue;
        }
    }
}

int RankTwoRelaxationMaxCut::getMaxCutValue() const {
    return bestCutValue;
}

const std::vector<int>& RankTwoRelaxationMaxCut::getBestSet() const {
    return bestSet;
}

}  // namespace Koala
