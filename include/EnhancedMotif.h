#ifndef ENHANCEDMOTIF_H
#define ENHANCEDMOTIF_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>
#include "MotifGraphBuilder.h"
#include "TabuSearchMotif.h"
#include "DNASequence.h"
#include "Graph.h"

// Helper structure for motif debugging
struct MotifDebugInfo {
    std::string motif;
    std::vector<int> positions;
    std::vector<std::string> sequences;
    int qualityThreshold;

    void print() const;
};

class EnhancedMotifGraphBuilder : public MotifGraphBuilder {
private:
    MotifDebugInfo debugInfo;
    bool isValidVertex(const KmerOccurrence &occ) const;
    bool shouldAddEdge(const KmerOccurrence &occ1, const KmerOccurrence &occ2) const;
    void verifyInjectedMotifs(const Graph &g) const;

public:
    EnhancedMotifGraphBuilder(
        const std::vector<DNASequence> &seqs,
        int kMin_,
        int kMax_,
        int qualityThreshold_,
        int positionMultiplier_,
        int allowedMismatch_,
        const MotifDebugInfo &debug);

    Graph build() const override;
};

class EnhancedTabuSearchMotif : public TabuSearchMotif {
public:
    using TabuSearchMotif::TabuSearchMotif;

protected:
    int evaluateMotifSolution(const Solution& S) const override;
    void diversify(Solution& S) override;
};

#endif //ENHANCEDMOTIF_H