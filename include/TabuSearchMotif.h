//
// Created by konrad_guest on 18/01/2025.
//
#ifndef TABUSEARCHMOTIF_H
#define TABUSEARCHMOTIF_H

#include "TabuSearchBase.h"

// A specialized TabuSearch that ensures "exactly one vertex per sequence."
// If we have S sequences, we want a clique with size = S.
class TabuSearchMotif : public TabuSearchBase {
private:
    int numSequences; // total # of sequences => we want a clique of this size

    void initialize() override;
    Solution selectBestNeighbor(const Solution &S, const std::vector<Solution> &neighbors) override;

    // Check if a solution has at most one vertex from each sequence
    bool validMotifSolution(const Solution &S) const;

public:
    TabuSearchMotif(const Graph &G, int T1_sz, int T2_sz, int maxIter, int numSeq);

    // We override run() to stop early once we get a solution of size = numSequences
    Solution run() override;
};

#endif // TABUSEARCHMOTIF_H
