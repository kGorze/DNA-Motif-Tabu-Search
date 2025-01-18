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
        throw std::runtime_error("Niepoprawny zakres k-merów (4 <= kMin <= kMax <= 9).");
    }
    if (positionMultiplier_ <= 0) {
        throw std::runtime_error("positionMultiplier musi być > 0.");
    }
}

std::vector<MotifGraphBuilder::KmerOccurrence> MotifGraphBuilder::collectKmers() const {
    std::vector<KmerOccurrence> occurrences;

    for (int seqIndex = 0; seqIndex < (int)sequences.size(); seqIndex++) {
        const DNASequence &seq = sequences[seqIndex];
        const std::string &bases = seq.bases;
        const std::vector<int> &qual = seq.qualityScores;

        int len = (int)bases.size();
        for (int k = kMin; k <= kMax; k++) {
            if (k > len) break;

            for (int startPos = 0; startPos + k <= len; startPos++) {
                bool aboveThreshold = true;
                for (int offset = 0; offset < k; offset++) {
                    if (qual[startPos + offset] < qualityThreshold) {
                        aboveThreshold = false;
                        break;
                    }
                }
                if (!aboveThreshold) continue;

                std::string kmer = bases.substr(startPos, k);

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
    std::vector<KmerOccurrence> occs = collectKmers();
    int n = (int)occs.size();

    Graph motifGraph(n);

    // Uzupełniamy vertexInfo
    for (int i = 0; i < n; i++) {
        VertexData vd;
        vd.sequenceIndex = occs[i].seqIndex;
        vd.position = occs[i].position;
        vd.kmer = occs[i].kmer;
        motifGraph.vertexInfo[i] = vd;
    }

    // Dodajemy krawędzie
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (Graph::canConnect(motifGraph.vertexInfo[i], motifGraph.vertexInfo[j], positionMultiplier)) {
                motifGraph.add_edge(i, j);
            }
        }
    }

    return motifGraph;
}
