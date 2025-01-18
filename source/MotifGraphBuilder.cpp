#include "../include/MotifGraphBuilder.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

MotifGraphBuilder::MotifGraphBuilder(const std::vector<DNASequence>& seqs,
                                     int kMin_,
                                     int kMax_,
                                     int qualityThreshold_,
                                     int positionMultiplier_,
                                     int allowedMismatch_)
    : sequences(seqs),
      kMin(kMin_),
      kMax(kMax_),
      qualityThreshold(qualityThreshold_),
      positionMultiplier(positionMultiplier_),
      allowedMismatches(allowedMismatch_)
{
    if (kMin < 4 || kMax > 100 || kMin > kMax) {
        throw std::runtime_error("Invalid k range (4..9 recommended).");
    }
    if (positionMultiplier <= 0) {
        throw std::runtime_error("positionMultiplier must be > 0!");
    }
    if (allowedMismatches < 0) {
        throw std::runtime_error("allowedMismatches cannot be negative!");
    }
}

static int countDifferences(const std::string& k1, const std::string& k2) {
    if (k1.size() != k2.size()) return 999999;
    int diff = 0;
    for (size_t i = 0; i < k1.size(); i++) {
        if (k1[i] != k2[i]) diff++;
        if (diff > 2) break;
    }
    return diff;
}

std::vector<MotifGraphBuilder::KmerOccurrence> 
MotifGraphBuilder::collectKmers() const {
    std::vector<KmerOccurrence> occs;
    int removedCount = 0;

    for (int seqIndex = 0; seqIndex < (int)sequences.size(); seqIndex++) {
        const DNASequence& seq = sequences[seqIndex];
        const std::string& bases = seq.bases;
        const std::vector<int>& quals = seq.qualityScores;
        int len = (int)bases.size();

        for (int k = kMin; k <= kMax; k++) {
            if (k > len) break;
            for (int startPos = 0; startPos + k <= len; startPos++) {
                bool pass = true;
                for (int offset = 0; offset < k; offset++) {
                    if (quals[startPos + offset] < qualityThreshold) {
                        pass = false;
                        removedCount++;
                        break;
                    }
                }
                if (!pass) continue;

                KmerOccurrence occ;
                occ.seqIndex = seqIndex;
                occ.position = startPos;
                occ.kmer = bases.substr(startPos, k);
                occs.push_back(occ);
            }
        }
    }

    return occs;
}

std::vector<MotifGraphBuilder::KmerOccurrence> 
MotifGraphBuilder::filterRepeatingPatterns(const std::vector<KmerOccurrence>& occs) const {
    return occs;
}

Graph MotifGraphBuilder::build() const {
    std::vector<KmerOccurrence> occs = collectKmers();
    auto filtered_occs = filterRepeatingPatterns(occs);

    int n = (int)filtered_occs.size();
    Graph motifGraph(n);
    motifGraph.vertexInfo.resize(n);

    for (int i = 0; i < n; i++) {
        motifGraph.vertexInfo[i].sequenceIndex = filtered_occs[i].seqIndex;
        motifGraph.vertexInfo[i].position = filtered_occs[i].position;
        motifGraph.vertexInfo[i].kmer = filtered_occs[i].kmer;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (filtered_occs[i].seqIndex != filtered_occs[j].seqIndex) {
                int diffs = countDifferences(filtered_occs[i].kmer, filtered_occs[j].kmer);
                if (diffs <= allowedMismatches) {
                    int posDiff = std::abs(
                        filtered_occs[i].position - 
                        filtered_occs[j].position
                    );
                    int lenK = (int)filtered_occs[i].kmer.size();
                    if (posDiff <= positionMultiplier * lenK) {
                        motifGraph.add_edge(i, j);
                    }
                }
            }
        }
    }

    return motifGraph;
}
