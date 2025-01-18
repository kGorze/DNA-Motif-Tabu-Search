#ifndef MOTIFGRAPHBUILDER_H
#define MOTIFGRAPHBUILDER_H

#include <vector>
#include "DNASequence.h"
#include "Graph.h"

// Builds a graph suitable for motif finding:
//   - Each vertex corresponds to a valid k-mer occurrence in one sequence
//   - Edges connect vertices if they share the same k-mer, are in different sequences,
//     and satisfy the position difference constraint
//
// Usage example:
//   MotifGraphBuilder builder(seqs, kMin, kMax, qualityThreshold, positionMultiplier);
//   Graph motifGraph = builder.build();

class MotifGraphBuilder {
private:
    const std::vector<DNASequence> &sequences;
    int kMin;
    int kMax;
    int qualityThreshold;
    int positionMultiplier; // e.g. 10 => position difference <= 10*k

    // A temporary struct for storing the (seqIndex, position, k-mer) before building the Graph
    struct KmerOccurrence {
        int seqIndex;
        int position;
        std::string kmer;
    };

    // Collect all valid k-mer occurrences across all sequences
    std::vector<KmerOccurrence> collectKmers() const;

public:
    MotifGraphBuilder(
        const std::vector<DNASequence> &seqs,
        int kMin_,
        int kMax_,
        int qualityThreshold_,
        int positionMultiplier_
    );

    Graph build() const;
};

#endif // MOTIFGRAPHBUILDER_H
