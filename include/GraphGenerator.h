//
// Created by konrad_guest on 16/12/2024.
//

#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include "Graph.h"

class GraphGenerator {
public:
    // Generate a graph using the p-generator:
    // Inputs: n (number of vertices), a,b (two real numbers with 0 < a < b < 1)
    // Procedure:
    // w[i] = uniform(a,b) for i=1 to n
    // For each pair (i,j), i<j, add edge with probability (w[i]+w[j])/2
    static Graph generate(int n, double a, double b);
};


#endif

