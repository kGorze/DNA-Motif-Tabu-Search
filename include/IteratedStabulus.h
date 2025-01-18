#ifndef ITERATEDSTABULUS_H
#define ITERATEDSTABULUS_H

#include "Graph.h"
#include "Stabulus.h"

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
