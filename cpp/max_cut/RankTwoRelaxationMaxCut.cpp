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

void RankTwoRelaxationMaxCut::distributeThetaEvenly() {
    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        theta[i] = 2 * M_PI * i / graph->numberOfNodes();
    }
}

std::vector<bool> RankTwoRelaxationMaxCut::procedureCut() {
    double bestValue = -std::numeric_limits<double>::infinity();
    std::vector<bool> bestCut(graph->numberOfNodes()), x(graph->numberOfNodes());

    for (double alpha = 0; alpha <= M_PI; alpha += 0.01) {
        for (int i = 0; i < graph->numberOfNodes(); ++i) {
            x[i] = (theta[i] >= alpha && theta[i] < alpha + M_PI) ? true : false;
        }

        double value = calculateCutValue(x);
        if (value > bestValue) {
            bestValue = value;
            bestCut = x;
        }
    }
    return bestCut;
}

std::vector<double> RankTwoRelaxationMaxCut::calculateGradient(const std::vector<double>& theta) {
    std::vector<double> gradient(graph->numberOfNodes(), 0.0);
    graph->forEdges([&](NetworKit::node j, NetworKit::node k, NetworKit::edgeweight w) {
        gradient[j] += w * sin(theta[k] - theta[j]);
        gradient[k] += w * sin(theta[j] - theta[k]);
    });
    return gradient;
}


void RankTwoRelaxationMaxCut::gradientDescent(double alpha, int maxIterations) {
    for (int iter = 0; iter < maxIterations; ++iter) {
        std::vector<double> gradient = calculateGradient(theta);
        for (int j = 0; j < graph->numberOfNodes(); ++j) {
            theta[j] -= alpha * gradient[j];
        }
    }

    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        if (theta[i] < 0) theta[i] += 2 * M_PI;
        if (theta[i] >= 2 * M_PI) theta[i] -= 2 * M_PI;
    }

    std::sort(theta.begin(), theta.end());
}

void RankTwoRelaxationMaxCut::perturbTheta() {
    distributeThetaEvenly();
    gradientDescent(alpha, maxIterations);
}

void RankTwoRelaxationMaxCut::run() {
    theta.resize(graph->numberOfNodes());
    perturbTheta();
    maxCutSet = procedureCut();
    maxCutValue = calculateCutValue(maxCutSet);

    if (maxCutSet[graph->numberOfNodes() - 1] == true) {
        maxCutSet[graph->numberOfNodes() - 1] = false;
        double candidateValue = calculateCutValue(maxCutSet);
        if (candidateValue > maxCutValue) {
            maxCutValue = candidateValue;
        }
    }
}

}  // namespace Koala
