//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/MotifGraphBuilder.h"
#include <stdexcept>

MotifGraphBuilder::MotifGraphBuilder(
    const std::vector<DNASequence> &seqs,
    int kMin_,
    int kMax_,
    int qualityThreshold_,
    int positionMultiplier_
)
: sequences(seqs),
  kMin(kMin_),
  kMax(kMax_),
  qualityThreshold(qualityThreshold_),
  positionMultiplier(positionMultiplier_)
{
    if (kMin < 4 || kMax > 9 || kMin > kMax) {
        throw std::runtime_error("Invalid k-mer length range. Must be between 4 and 9.");
    }
    if (positionMultiplier_ <= 0) {
        throw std::runtime_error("positionMultiplier must be positive.");
    }
}

std::vector<MotifGraphBuilder::KmerOccurrence> MotifGraphBuilder::collectKmers() const {
    std::vector<KmerOccurrence> occurrences;

    for (int seqIndex = 0; seqIndex < (int)sequences.size(); seqIndex++) {
        const DNASequence &seq = sequences[seqIndex];
        const std::string &bases = seq.bases;
        const std::vector<int> &qual = seq.qualityScores;

        int len = (int)bases.size();
        // For each possible k in [kMin, kMax]
        for (int k = kMin; k <= kMax; k++) {
            if (k > len) break;

            for (int startPos = 0; startPos + k <= len; startPos++) {
                // Check if all positions in this k-mer meet the quality threshold
                bool aboveThreshold = true;
                for (int offset = 0; offset < k; offset++) {
                    if (qual[startPos + offset] < qualityThreshold) {
                        aboveThreshold = false;
                        break;
                    }
                }
                if (!aboveThreshold) continue;

                // Extract the k-mer
                std::string kmer = bases.substr(startPos, k);

                // Store occurrence
                KmerOccurrence occ;
                occ.seqIndex = seqIndex;
                occ.position = startPos;
                occ.kmer = kmer;
                occurrences.push_back(occ);
            }
        }
    }

    return occurrences;
}

Graph MotifGraphBuilder::build() const {
    // 1) Collect all valid k-mer occurrences
    std::vector<KmerOccurrence> occs = collectKmers();
    int n = (int)occs.size();

    // 2) Create a Graph with n vertices
    Graph motifGraph(n);

    // 3) Fill vertexInfo
    for (int i = 0; i < n; i++) {
        VertexData vd;
        vd.sequenceIndex = occs[i].seqIndex;
        vd.position = occs[i].position;
        vd.kmer = occs[i].kmer;
        motifGraph.vertexInfo[i] = vd;
    }

    // 4) Add edges according to motif rule
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // check canConnect
            if (Graph::canConnect(motifGraph.vertexInfo[i], motifGraph.vertexInfo[j], positionMultiplier)) {
                motifGraph.add_edge(i, j);
            }
        }
    }

    return motifGraph;
}
