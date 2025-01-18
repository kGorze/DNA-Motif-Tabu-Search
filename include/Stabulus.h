#ifndef STABULUS_H
#define STABULUS_H

#include "Graph.h"
#include <vector>

// Prosty algorytm do znajdowania stabilnego zbioru wielkości k
class Stabulus {
private:
    const Graph &G;
    int k;

public:
    Stabulus(const Graph &graph, int k_);
    bool findStableSetOfSizeK(std::vector<int> &stableSet);
};

#endif
