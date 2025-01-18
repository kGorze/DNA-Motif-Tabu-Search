#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include "Graph.h"

class GraphGenerator {
public:
    // Generuje graf używając w[i] = uniform(a,b) i prawdopodobieństwa (w[i]+w[j])/2
    static Graph generate(int n, double a, double b);
};

#endif
