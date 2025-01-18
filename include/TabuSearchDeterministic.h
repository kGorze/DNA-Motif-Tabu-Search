#ifndef TABUSEARCHDETERMINISTIC_H
#define TABUSEARCHDETERMINISTIC_H

#include "TabuSearchBase.h"

class TabuSearchDeterministic : public TabuSearchBase {
private:
    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, 
                                const std::vector<Solution> &neighbors) override;

public:
    TabuSearchDeterministic(const Graph &G, int T1_sz, int T2_sz, int maxIter);
    Solution run() override;
};

#endif // TABUSEARCHDETERMINISTIC_H
