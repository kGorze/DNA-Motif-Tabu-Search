#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include "Graph.h"

class GraphGenerator {
public:
    // Generuje graf używając metody p-generator:
    // w[i] = uniform(a,b)
    // Krawędź (i,j) powstaje z prawdopodobieństwem (w[i]+w[j])/2
    static Graph generate(int n, double a, double b);
};

#endif
