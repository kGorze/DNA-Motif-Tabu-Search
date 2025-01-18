#include "../include/Stabulus.h"

Stabulus::Stabulus(const Graph &graph, int k_): G(graph), k(k_) {}

bool Stabulus::findStableSetOfSizeK(std::vector<int> &stableSet) {
    stableSet.clear();
    for (int v = 0; v < G.n; v++) {
        bool canAdd = true;
        for (int u : stableSet) {
            if (G.is_edge(u,v)) {
                canAdd = false; 
                break;
            }
        }
        if (canAdd) {
            stableSet.push_back(v);
            if ((int)stableSet.size() == k) return true;
        }
    }
    return ((int)stableSet.size() == k);
}
