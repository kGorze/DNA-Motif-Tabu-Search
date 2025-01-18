#include "../include/IteratedStabulus.h"

IteratedStabulus::IteratedStabulus(const Graph &graph, int initialK_, int maxFail_)
    : G(graph), initialK(initialK_), maxFail(maxFail_) {}

int IteratedStabulus::runFindMaxStableSet(std::vector<int> &bestSet) {
    int k = initialK;
    int failCount = 0;
    int bestK = initialK;
    bestSet.clear();

    while (failCount < maxFail) {
        Stabulus stab(G, k);
        std::vector<int> candidate;
        bool found = stab.findStableSetOfSizeK(candidate);
        if (found) {
            bestSet = candidate;
            bestK = k;
            k++;
            failCount = 0;
        } else {
            failCount++;
            k++;
        }
    }

    return bestK;
}
