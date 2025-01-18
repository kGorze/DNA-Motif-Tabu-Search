//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/EnhancedMotif.h"

// Global RNG for the tabu search
static std::mt19937 g_rng2((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

static int countDifferences(const std::string &s1, const std::string &s2) {
    if (s1.size() != s2.size()) return 999999;
    int diff = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) diff++;
    }
    return diff;
}

// MotifDebugInfo implementation
void MotifDebugInfo::print() const {
    std::cout << "\n=== Motif Debug Info ===\n"
              << "Injected motif: " << motif << "\n"
              << "Positions: ";
    for (int pos : positions) {
        std::cout << pos << " ";
    }
    std::cout << "\nSequences at injection points:\n";
    for (const auto &seq : sequences) {
        std::cout << "  " << seq << "\n";
    }
}

// EnhancedMotifGraphBuilder implementations
EnhancedMotifGraphBuilder::EnhancedMotifGraphBuilder(
    const std::vector<DNASequence> &seqs,
    int kMin_,
    int kMax_,
    int qualityThreshold_,
    int positionMultiplier_,
    int allowedMismatch_,
    const MotifDebugInfo &debug)
    : MotifGraphBuilder(seqs, kMin_, kMax_, qualityThreshold_, 
                       positionMultiplier_, allowedMismatch_),
      debugInfo(debug) {}

bool EnhancedMotifGraphBuilder::isValidVertex(const KmerOccurrence &occ) const {
    const DNASequence &seq = sequences[occ.seqIndex];
    for (int i = 0; i < (int)occ.kmer.length(); i++) {
        if (seq.qualityScores[occ.position + i] < qualityThreshold) {
            return false;
        }
    }
    if (occ.position < MIN_START_POS) {
        return false;
    }
    return true;
}

bool EnhancedMotifGraphBuilder::shouldAddEdge(
    const KmerOccurrence &occ1, 
    const KmerOccurrence &occ2) const {
    if (occ1.seqIndex == occ2.seqIndex) {
        return false;
    }
    int diffs = countDifferences(occ1.kmer, occ2.kmer);
    if (diffs > allowedMismatches) {
        return false;
    }
    int posDiff = std::abs(occ1.position - occ2.position);
    int maxAllowedDiff = positionMultiplier * (int)occ1.kmer.length();
    if (posDiff > maxAllowedDiff) {
        return false;
    }
    return true;
}

void EnhancedMotifGraphBuilder::verifyInjectedMotifs(const Graph &g) const {
    std::cout << "\n=== Verifying Injected Motifs ===\n";

    std::vector<int> motifVertices;
    for (int i = 0; i < g.n; i++) {
        const auto &vInfo = g.vertexInfo[i];
        for (size_t j = 0; j < debugInfo.positions.size(); j++) {
            if (vInfo.sequenceIndex == (int)j &&
                vInfo.position == debugInfo.positions[j] &&
                countDifferences(vInfo.kmer, debugInfo.motif) <= allowedMismatches) {
                motifVertices.push_back(i);
                std::cout << "Found injected motif vertex: seq=" << j 
                          << " pos=" << vInfo.position 
                          << " kmer=" << vInfo.kmer << "\n";
            }
        }
    }

    std::cout << "\nChecking edges between motif vertices:\n";
    for (size_t i = 0; i < motifVertices.size(); i++) {
        for (size_t j = i + 1; j < motifVertices.size(); j++) {
            int v1 = motifVertices[i];
            int v2 = motifVertices[j];
            if (g.is_edge(v1, v2)) {
                std::cout << "Edge exists between motifs: " 
                          << v1 << " - " << v2 << "\n";
            } else {
                std::cout << "Missing edge between motifs: " 
                          << v1 << " - " << v2 << "\n";
            }
        }
    }
}

Graph EnhancedMotifGraphBuilder::build() const {
    auto kmers = collectKmers();
    std::cout << "Collected " << kmers.size() << " raw k-mers\n";

    std::vector<KmerOccurrence> filtered;
    filtered.reserve(kmers.size());
    for (auto &occ : kmers) {
        if (isValidVertex(occ)) {
            filtered.push_back(occ);
        }
    }
    std::cout << "After isValidVertex check: " << filtered.size() << " k-mers\n";

    auto filtered2 = filterRepeatingPatterns(filtered);
    std::cout << "After filterRepeatingPatterns: " << filtered2.size() << " k-mers\n";

    Graph g((int)filtered2.size());
    g.vertexInfo.resize((int)filtered2.size());

    for (int i = 0; i < (int)filtered2.size(); i++) {
        g.vertexInfo[i].sequenceIndex = filtered2[i].seqIndex;
        g.vertexInfo[i].position = filtered2[i].position;
        g.vertexInfo[i].kmer = filtered2[i].kmer;
    }

    int edgeCount = 0;
    for (int i = 0; i < (int)filtered2.size() - 1; i++) {
        for (int j = i + 1; j < (int)filtered2.size(); j++) {
            if (shouldAddEdge(filtered2[i], filtered2[j])) {
                g.add_edge(i, j);
                edgeCount++;
            }
        }
    }
    
    std::cout << "Built graph with " << g.n << " vertices and "
              << edgeCount << " edges\n";

    verifyInjectedMotifs(g);

    return g;
}

// EnhancedTabuSearchMotif implementations
int EnhancedTabuSearchMotif::evaluateMotifSolution(const Solution &S) const {
    int score = S.size * 1000;

    std::vector<int> vertices;
    vertices.reserve(S.size);
    for (int i = 0; i < G.n; i++) {
        if (S.inClique[i]) {
            vertices.push_back(i);
        }
    }
    
    int totalDiff = 0;
    int pairs = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i + 1; j < vertices.size(); j++) {
            const auto &k1 = G.vertexInfo[vertices[i]].kmer;
            const auto &k2 = G.vertexInfo[vertices[j]].kmer;
            totalDiff += countDifferences(k1, k2);
            pairs++;
        }
    }
    if (pairs > 0) {
        double avgDiff = (double)totalDiff / pairs;
        score -= (int)(avgDiff * 20.0);
    }

    int totalPosDiff = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i + 1; j < vertices.size(); j++) {
            int pos1 = G.vertexInfo[vertices[i]].position;
            int pos2 = G.vertexInfo[vertices[j]].position;
            totalPosDiff += std::abs(pos1 - pos2);
        }
    }
    if (pairs > 0) {
        double avgPosDiff = (double)totalPosDiff / pairs;
        score -= (int)(avgPosDiff * 5.0);
    }
    
    return score;
}

void EnhancedTabuSearchMotif::diversify(Solution &S) {
    if (itersSinceImprovement > MaxIter / 2) {
        std::vector<int> inClique;
        for (int i = 0; i < G.n; i++) {
            if (S.inClique[i]) {
                inClique.push_back(i);
            }
        }
        if (!inClique.empty()) {
            std::shuffle(inClique.begin(), inClique.end(), g_rng2);
            int toRemove = std::min(3, (int)inClique.size());
            for (int i = 0; i < toRemove; i++) {
                S.inClique[inClique[i]] = false;
            }
            S.size -= toRemove;
        }
    }
}