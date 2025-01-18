//
// Created by konrad_guest on 16/12/2024.
//

#ifndef TABUSEARCHDETERMINISTIC_H
#define TABUSEARCHDETERMINISTIC_H

#include "TabuSearchBase.h"

class TabuSearchDeterministic : public TabuSearchBase {
private:
    // Deterministic variant initialization and move selection
    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, const std::vector<Solution> &neighbors) override;

public:
    TabuSearchDeterministic(const Graph &G, int T1_sz, int T2_sz, int maxIter);
    Solution run() override;
};

#endif
