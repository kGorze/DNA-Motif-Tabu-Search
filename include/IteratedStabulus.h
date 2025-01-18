//
// Created by konrad_guest on 16/12/2024.
//
#ifndef ITERATEDSTABULUS_H
#define ITERATEDSTABULUS_H

#include "Graph.h"
#include "Stabulus.h"

// Iterated-Stabulus tries k, k+1, k+2,... until fail or time limit
class IteratedStabulus {
private:
    const Graph &G;
    int initialK;
    int maxFail;

public:
    IteratedStabulus(const Graph &graph, int initialK_, int maxFail_);
    int runFindMaxStableSet(std::vector<int> &bestSet);
};

#endif
