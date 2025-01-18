//
// Created by konrad_guest on 16/12/2024.
//

#ifndef TABUSEARCHPROBABILISTIC_H
#define TABUSEARCHPROBABILISTIC_H

#include "TabuSearchBase.h"

class TabuSearchProbabilistic : public TabuSearchBase {
private:
    int k; // number of vertices to remove in shake-up
    int L; // sample size

    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, const std::vector<Solution> &neighbors) override;
    void randomShakeUp(Solution &S, int k);
    std::vector<Solution> sampleAugmentingNeighbors(const std::vector<Solution> &N_plus, int L);

public:
    TabuSearchProbabilistic(const Graph &G, int T1_sz, int T2_sz, int maxIter, int k_, int L_);
    Solution run() override;
};

#endif
