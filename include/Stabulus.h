//
// Created by konrad_guest on 16/12/2024.
//

#ifndef STABULUS_H
#define STABULUS_H

#include "Graph.h"
#include <vector>

// STABULUS: searches a stable set of size k in graph G
// Not strictly optimizing, just tries to find a stable set of given size k.

class Stabulus {
private:
    const Graph &G;
    int k; // the size of stable set to find

public:
    Stabulus(const Graph &graph, int k_);
    bool findStableSetOfSizeK(std::vector<int> &stableSet);
};

#endif
