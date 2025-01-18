#ifndef TABUSEARCHPROBABILISTIC_H
#define TABUSEARCHPROBABILISTIC_H

#include "TabuSearchBase.h"

class TabuSearchProbabilistic : public TabuSearchBase {
private:
    int k;  // liczba wierzcholkow usuwanych w random shake-up
    int L;  // rozmiar probkowania N+(S)

    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, 
                                const std::vector<Solution> &neighbors) override;
    void randomShakeUp(Solution &S, int kk);
    std::vector<Solution> sampleAugmentingNeighbors(const std::vector<Solution> &N_plus, int LL);

public:
    TabuSearchProbabilistic(const Graph &G, int T1_sz, int T2_sz, int maxIter, int k_, int L_);
    Solution run() override;
};

#endif // TABUSEARCHPROBABILISTIC_H
