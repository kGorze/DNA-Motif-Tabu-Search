#ifndef TABUSEARCHMOTIF_H
#define TABUSEARCHMOTIF_H

#include "TabuSearchBase.h"

// Specjalizowana wersja Tabu Search do znajdowania kliki,
// w której jest dokładnie jeden wierzchołek z każdej sekwencji.
// Jeśli mamy S sekwencji, chcemy kliki o rozmiarze = S.
class TabuSearchMotif : public TabuSearchBase {
private:
    int numSequences; // liczba sekwencji

    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, const std::vector<Solution> &neighbors) override;
    bool validMotifSolution(const Solution &S) const;

public:
    TabuSearchMotif(const Graph &G, int T1_sz, int T2_sz, int maxIter, int numSeq);
    Solution run() override;
};

#endif // TABUSEARCHMOTIF_H
