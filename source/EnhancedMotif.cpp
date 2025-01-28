//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/EnhancedMotif.h"
#include <map>
#include <set>

// Global RNG for the tabu search
static std::mt19937 g_rng2((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// Add this structure at the top of the file, after the includes
struct MotifData {
    std::vector<std::pair<int, int>> positions;
    bool isInjected;
    bool isRelatedToInjected;
    
    MotifData() : isInjected(false), isRelatedToInjected(false) {}
    
    // Add comparison operator
    bool operator<(const MotifData& other) const {
        if (positions.size() != other.positions.size())
            return positions.size() < other.positions.size();
        if (isInjected != other.isInjected)
            return isInjected < other.isInjected;
        if (isRelatedToInjected != other.isRelatedToInjected)
            return isRelatedToInjected < other.isRelatedToInjected;
        return positions < other.positions;
    }
};

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
    std::cout << "\n=== Verifying Injected Motifs ===\n\n";
    
    // First show any found cliques/motifs
    std::cout << "Found k-mer patterns:\n";
    std::map<std::string, std::vector<int>> kmerGroups;  // kmer -> list of vertex indices
    
    for (int v = 0; v < g.n; v++) {
        const auto& info = g.vertexInfo[v];
        kmerGroups[info.kmer].push_back(v);
    }
    
    // Show k-mers that appear in multiple sequences
    for (const auto& group : kmerGroups) {
        std::set<int> seqIndices;
        for (int v : group.second) {
            seqIndices.insert(g.vertexInfo[v].sequenceIndex);
        }
        
        if (seqIndices.size() >= 3) {  // Show if present in at least 3 sequences
            std::cout << "\nMotif pattern: " << group.first << "\n";
            std::cout << "Found in " << seqIndices.size() << " sequences:\n";
            for (int v : group.second) {
                const auto& info = g.vertexInfo[v];
                std::cout << "  Seq " << info.sequenceIndex + 1 
                         << " at pos " << info.position << "\n";
            }
        }
    }
    
    // If we have an injected motif, verify it as well
    if (!debugInfo.motif.empty()) {
        std::cout << "\nChecking for injected motif: " << debugInfo.motif << "\n";
        std::vector<int> motifVertices;
        for (int v = 0; v < g.n; v++) {
            const auto& info = g.vertexInfo[v];
            if (info.position == debugInfo.positions[info.sequenceIndex]) {
                motifVertices.push_back(v);
            }
        }
        
        if (!motifVertices.empty()) {
            std::cout << "Found " << motifVertices.size() 
                      << " vertices matching injected positions\n";
        }
    }
    
    std::cout << "\nVerification complete.\n";
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

void verifyInjectedMotif(const std::vector<DNASequence>& sequences, const std::string& injectedMotif) {
    std::cout << "\n=== Verifying Injected Motifs ===\n\n";
    
    // Group k-mers by length first
    std::map<int, std::map<std::string, MotifData>> kmersByLength;
    
    // First pass: collect all kmers and group by length
    for (int len = 4; len <= 9; len++) {
        for (int seqIdx = 0; seqIdx < sequences.size(); seqIdx++) {
            const std::string& seq = sequences[seqIdx].bases;
            for (int i = 0; i <= seq.length() - len; i++) {
                std::string kmer = seq.substr(i, len);
                kmersByLength[len][kmer].positions.push_back({seqIdx + 1, i + 1});
            }
        }
    }
    
    // Process each length group
    for (const auto& lengthGroup : kmersByLength) {
        int currentLength = lengthGroup.first;
        std::cout << "\n=== Length " << currentLength << " Motifs ===\n\n";
        
        // Sort kmers within this length group
        std::vector<std::pair<std::string, MotifData>> sortedKmers = {lengthGroup.second.begin(), lengthGroup.second.end()};
        std::sort(sortedKmers.begin(), sortedKmers.end());
        
        // Merge duplicate kmers
        std::map<std::string, MotifData> mergedKmers;
        for (const auto& kmerData : sortedKmers) {
            mergedKmers[kmerData.first] = kmerData.second;
        }
        
        // Output sorted, merged kmers
        for (const auto& kmerData : mergedKmers) {
            const std::string& kmer = kmerData.first;
            const MotifData& data = kmerData.second;
            
            // Count unique sequences
            std::set<int> uniqueSeqs;
            for (const auto& pos : data.positions) {
                uniqueSeqs.insert(pos.first);
            }
            
            if (uniqueSeqs.size() >= 3) {  // Only show kmers present in at least 3 sequences
                if (kmer == injectedMotif) {
                    std::cout << "Motif pattern: " << kmer << " (injected motif)\n";
                    mergedKmers[kmer].isInjected = true;
                } else if (kmer.find(injectedMotif) != std::string::npos || 
                          injectedMotif.find(kmer) != std::string::npos) {
                    std::cout << "Motif pattern: " << kmer << " (related to injected motif)\n";
                    mergedKmers[kmer].isRelatedToInjected = true;
                } else {
                    std::cout << "Motif pattern: " << kmer << " (natural motif)\n";
                }
                std::cout << "Found in " << uniqueSeqs.size() << " sequences:\n";
                for (const auto& pos : data.positions) {
                    std::cout << "  Seq " << pos.first << " at pos " << pos.second << "\n";
                }
                std::cout << "\n";
            }
        }
    }
    
    // Add this section at the end before "Verification complete"
    std::cout << "\n=== SUMMARY OF INJECTED MOTIF PATTERNS ===\n";
    for (const auto& lengthGroup : kmersByLength) {
        for (const auto& kmerData : lengthGroup.second) {
            const std::string& kmer = kmerData.first;
            const MotifData& data = kmerData.second;
            
            if (data.isInjected || data.isRelatedToInjected) {
                std::cout << "\nMotif pattern: " << kmer 
                         << (data.isInjected ? " (INJECTED)" : " (related to injected)")
                         << "\nFound in " << data.positions.size() << " sequences:\n";
                for (const auto& pos : data.positions) {
                    std::cout << "  Seq " << pos.first << " at pos " << pos.second << "\n";
                }
            }
        }
    }
    std::cout << "\n=== End of Injected Motif Summary ===\n\n";
    
    // Add this section before "Verification complete"
    std::cout << "\n=== MOST COMMON MOTIFS ===\n";
    
    // Store all motifs with their occurrence count and positions
    std::map<std::string, std::map<int, int>> motifOccurrences; // motif -> {seqIdx -> pos}
    
    // Collect all 9-mers
    for (const auto& seq : sequences) {
        for (int i = 0; i <= seq.bases.length() - 9; i++) {
            std::string kmer = seq.bases.substr(i, 9);
            motifOccurrences[kmer][&seq - &sequences[0] + 1] = i + 1;
        }
    }
    
    // Sort motifs by occurrence count
    std::vector<std::pair<std::string, std::map<int, int>>> sortedMotifs;
    for (const auto& motif : motifOccurrences) {
        if (motif.second.size() >= 3) { // Only include motifs present in at least 3 sequences
            sortedMotifs.push_back(motif);
        }
    }
    
    std::sort(sortedMotifs.begin(), sortedMotifs.end(),
        [](const auto& a, const auto& b) {
            if (a.second.size() != b.second.size())
                return a.second.size() > b.second.size();
            return a.first < b.first;
        });
    
    // Show top motifs
    int motifCount = 0;
    for (const auto& motif : sortedMotifs) {
        std::cout << "\nMotif pattern: " << motif.first << "\n";
        std::cout << "Found in " << motif.second.size() << " sequences:\n";
        for (const auto& pos : motif.second) {
            std::cout << "  Seq " << pos.first << " at pos " << pos.second << "\n";
        }
        motifCount++;
    }
    
    std::cout << "\nTotal number of motifs found: " << motifCount << "\n";
    std::cout << "\nVerification complete.\n";
}