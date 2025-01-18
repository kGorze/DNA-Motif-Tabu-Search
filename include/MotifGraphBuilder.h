#ifndef MOTIFGRAPHBUILDER_H
#define MOTIFGRAPHBUILDER_H

#include <vector>
#include "DNASequence.h"
#include "Graph.h"

// Budowanie grafu do motif-finding:
// - Każdy wierzchołek to wystąpienie k-mera (o długości w zakresie [kMin, kMax])
//   w jednej z sekwencji, z odpowiednimi jakościami
// - Krawędzie łączą wierzchołki, jeśli k-mery są identyczne,
//   pochodzą z różnych sekwencji, a różnica pozycji ≤ positionMultiplier * k
class MotifGraphBuilder {
private:
    const std::vector<DNASequence> &sequences;
    int kMin;
    int kMax;
    int qualityThreshold;
    int positionMultiplier; // np. 10 => |pos1-pos2| ≤ 10*k

    struct KmerOccurrence {
        int seqIndex;
        int position;
        std::string kmer;
    };

    // Zbiera wszystkie k-mery spełniające próg jakości
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
