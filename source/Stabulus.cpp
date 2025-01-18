//
// Created by konrad_guest on 16/12/2024.
//

#include "../include/Stabulus.h"

// Dummy implementation: tries some heuristic to find stable set of size k
// This is just a placeholder for demonstration.
// A stable set in G is a clique in complement(G).
bool Stabulus::findStableSetOfSizeK(std::vector<int> &stableSet) {
    // Very naive approach: try a simple heuristic
    // Start empty and add vertices that don't conflict with chosen ones
    stableSet.clear();
    for (int v=0; v<G.n; v++) {
        bool canAdd = true;
        for (int u : stableSet) {
            if (G.adj[u][v]) { // edge means not independent
                canAdd=false;break;
            }
        }
        if(canAdd) {
            stableSet.push_back(v);
            if((int)stableSet.size()==k) return true;
        }
    }
    return (int)stableSet.size()==k;
}

Stabulus::Stabulus(const Graph &graph, int k_):G(graph),k(k_){}
