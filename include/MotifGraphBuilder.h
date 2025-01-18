#ifndef MOTIFGRAPHBUILDER_H
#define MOTIFGRAPHBUILDER_H

#include <vector>
#include <limits>
#include "DNASequence.h"
#include "Graph.h"

class MotifGraphBuilder {
protected:
    const std::vector<DNASequence>& sequences;
    int kMin;
    int kMax;
    int qualityThreshold;
    int positionMultiplier;
    int allowedMismatches;

    // Internal structure to represent k-mer occurrences
    struct KmerOccurrence {
        int seqIndex;
        int position;
        std::string kmer;
    };

    // Helper methods
    std::vector<KmerOccurrence> collectKmers() const;
    std::vector<KmerOccurrence> filterRepeatingPatterns(const std::vector<KmerOccurrence>& occs) const;
    
    static const int MIN_START_POS = 5; // Minimum allowed starting position

    int countDifferences(const std::string& kmer1, const std::string& kmer2) const {
        if (kmer1.length() != kmer2.length()) {
            return std::numeric_limits<int>::max(); // Różne długości - zwróć maksymalną wartość
        }

        int differences = 0;
        for (size_t i = 0; i < kmer1.length(); ++i) {
            // Porównuj znaki ignorując wielkość liter
            char c1 = std::toupper(kmer1[i]);
            char c2 = std::toupper(kmer2[i]);
            
            // Sprawdź czy nukleotydy są różne
            if (c1 != c2) {
                differences++;
                
                // Opcjonalna optymalizacja: jeśli przekroczyliśmy dozwoloną liczbę różnic,
                // możemy przerwać wcześniej
                if (differences > allowedMismatches) {
                    return differences;
                }
            }
        }
        
        return differences;
    }

public:
    // Constructor with parameters for building the motif graph
    MotifGraphBuilder(const std::vector<DNASequence>& seqs,
                     int kMin_,
                     int kMax_,
                     int qualityThreshold_,
                     int positionMultiplier_,
                     int allowedMismatch_ = 0);

    // Main method to construct the graph
    virtual Graph build() const;
};

#endif // MOTIFGRAPHBUILDER_H
