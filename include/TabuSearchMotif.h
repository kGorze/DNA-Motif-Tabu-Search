#ifndef TABUSEARCHMOTIF_H
#define TABUSEARCHMOTIF_H

#include "TabuSearchBase.h"
#include "Graph.h"

class TabuSearchMotif : public TabuSearchBase {
protected:
    int numSequences;

    // Nowe parametry do dywersyfikacji
    int diversificationInterval; 
    int diversificationCounter;

    // Override base class virtual methods
    void initialize() override;
    Solution selectBestNeighbor(const Solution& S,
                                const std::vector<Solution>& neighbors) override;

    // Helper methods
    bool validMotifSolution(const Solution& S) const;

    // Dodajemy metodę do lekkiej dywersyfikacji
    virtual int evaluateMotifSolution(const Solution& S) const;
    virtual void diversify(Solution& S);

public:
    TabuSearchMotif(const Graph& G, int T1_sz, int T2_sz, int maxIter, int numSeq);
    Solution run() override;
};

#endif // TABUSEARCHMOTIF_H
