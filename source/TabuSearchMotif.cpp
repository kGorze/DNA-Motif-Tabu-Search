#include "../include/TabuSearchMotif.h"
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>

static std::mt19937 g_rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

static int findAddedVertex(const Solution& oldSolution, const Solution& newSolution) {
    if (oldSolution.inClique.size() != newSolution.inClique.size()) {
        return -1;
    }
    int addedVertex = -1;
    int differences = 0;

    for (size_t i = 0; i < oldSolution.inClique.size(); ++i) {
        if (oldSolution.inClique[i] != newSolution.inClique[i]) {
            if (!oldSolution.inClique[i] && newSolution.inClique[i]) {
                addedVertex = (int)i;
            }
            differences++;
        }
    }
    return (differences == 1 && addedVertex != -1) ? addedVertex : -1;
}

static int countDifferences(const std::string &s1, const std::string &s2) {
    if (s1.size() != s2.size()) return 999999;
    int diff = 0;
    for (size_t i = 0; i < s1.size(); i++) {
        if (s1[i] != s2[i]) diff++;
    }
    return diff;
}

static Solution naiveMotifInit(const Graph &G, int numSeq) {
    Solution S;
    S.inClique.resize(G.n,false);
    S.size=0;

    // Try to find a densely connected region to start from
    std::vector<int> vertexDegrees(G.n);
    for (int i = 0; i < G.n; i++) {
        vertexDegrees[i] = 0;
        for (int j = 0; j < G.n; j++) {
            if (G.is_edge(i, j)) vertexDegrees[i]++;
        }
    }

    std::vector<bool> used(numSeq, false);
    std::vector<std::pair<int,int>> degreeVertices;
    for (int i = 0; i < G.n; i++) {
        degreeVertices.push_back({vertexDegrees[i], i});
    }
    std::sort(degreeVertices.rbegin(), degreeVertices.rend());

    // Try to start from high-degree vertices
    for (const auto &p : degreeVertices) {
        int v = p.second;
        int seqIdx = G.vertexInfo[v].sequenceIndex;
        if (!used[seqIdx]) {
            bool canPick = true;
            for (int u = 0; u < G.n; u++) {
                if (S.inClique[u]) {
                    if (!G.is_edge(u, v)) {
                        canPick=false;
                        break;
                    }
                }
            }
            if (canPick) {
                S.inClique[v] = true;
                S.size++;
                used[seqIdx] = true;
                if (S.size == numSeq) break;
            }
        }
    }
    return S;
}

TabuSearchMotif::TabuSearchMotif(const Graph& G, int T1_sz, int T2_sz, int maxIter, int numSeq)
    : TabuSearchBase(G, T1_sz, T2_sz, maxIter), 
      numSequences(numSeq),
      diversificationInterval(30),  // Reduced interval for more frequent diversification
      diversificationCounter(0)
{}

void TabuSearchMotif::initialize() {
    bestSol = naiveMotifInit(G, numSequences);
    bestFoundSize = bestSol.size;
    itersSinceImprovement = 0;

    T1_list.clear();
    T1_set.clear();
    T2_list.clear();
    T2_set.clear();

    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

bool TabuSearchMotif::validMotifSolution(const Solution& S) const {
    std::vector<int> used(numSequences, 0);
    std::vector<int> verts;

    for (int i = 0; i < G.n; i++) {
        if (S.inClique[i]) {
            int seqIdx = G.vertexInfo[i].sequenceIndex;
            used[seqIdx]++;
            if (used[seqIdx] > 1) return false;
            verts.push_back(i);
        }
    }
    
    // Check connectivity with some tolerance for temporary disconnections
    int allowedMissingEdges = S.size <= 3 ? 0 : 1;  // Allow 1 missing edge for larger cliques
    int missingEdges = 0;
    
    for (size_t i = 0; i < verts.size(); i++) {
        for (size_t j = i+1; j < verts.size(); j++) {
            if (!G.is_edge(verts[i], verts[j])) {
                missingEdges++;
                if (missingEdges > allowedMissingEdges) {
                    return false;
                }
            }
        }
    }
    return true;
}

int TabuSearchMotif::evaluateMotifSolution(const Solution& S) const {
    // Heavily weight size to prioritize finding larger cliques
    int score = S.size * 5000;

    std::vector<int> vertices;
    vertices.reserve(S.size);
    for (int i = 0; i < G.n; i++) {
        if (S.inClique[i]) {
            vertices.push_back(i);
        }
    }

    // Calculate sequence similarity score
    int totalDiff = 0;
    int pairs = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i+1; j < vertices.size(); j++) {
            const auto &k1 = G.vertexInfo[vertices[i]].kmer;
            const auto &k2 = G.vertexInfo[vertices[j]].kmer;
            totalDiff += countDifferences(k1, k2);
            pairs++;
        }
    }
    
    if (pairs > 0) {
        double avgDiff = (double)totalDiff / pairs;
        // Reduced penalty weight for differences
        score -= (int)(avgDiff * 5.0);
    }

    // Add bonus for position closeness
    int totalPosDiff = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i+1; j < vertices.size(); j++) {
            int pos1 = G.vertexInfo[vertices[i]].position;
            int pos2 = G.vertexInfo[vertices[j]].position;
            totalPosDiff += std::abs(pos1 - pos2);
        }
    }
    
    if (pairs > 0) {
        double avgPosDiff = (double)totalPosDiff / pairs;
        score -= (int)(avgPosDiff * 2.0);  // Small penalty for position differences
    }

    // Bonus for dense connectivity
    int edges = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i+1; j < vertices.size(); j++) {
            if (G.is_edge(vertices[i], vertices[j])) {
                edges++;
            }
        }
    }
    if (pairs > 0) {
        double density = (double)edges / pairs;
        score += (int)(density * 1000.0);  // Bonus for dense connections
    }

    return score;
}

Solution TabuSearchMotif::selectBestNeighbor(const Solution& S,
                                             const std::vector<Solution>& neighbors)
{
    int bestVal = -1;
    Solution bestNeighbor;
    bestNeighbor.size = -1;

    for (const auto& Sprime : neighbors) {
        if (!validMotifSolution(Sprime)) {
            continue;
        }
        if (T1_set.find(Sprime) != T1_set.end()) {
            continue;
        }
        bool augmenting = (Sprime.size > S.size);
        if (augmenting) {
            int addedV = findAddedVertex(S, Sprime);
            if (addedV != -1 && T2_set.find(addedV) != T2_set.end()) {
                if (Sprime.size <= bestFoundSize) {
                    continue;
                }
            }
        }

        int val = evaluateMotifSolution(Sprime);
        if (val > bestVal) {
            bestVal = val;
            bestNeighbor = Sprime;
        }
    }
    return bestNeighbor;
}

void TabuSearchMotif::diversify(Solution &S) {
    std::vector<int> inClique;
    for (int i = 0; i < G.n; i++) {
        if (S.inClique[i]) {
            inClique.push_back(i);
        }
    }
    if (inClique.empty()) return;

    std::shuffle(inClique.begin(), inClique.end(), g_rng);
    
    // Adaptive removal based on current solution size and search progress
    int toRemove;
    if (itersSinceImprovement < MaxIter / 3) {
        toRemove = 1;  // Light diversification early in search
    } else if (itersSinceImprovement < MaxIter * 2/3) {
        toRemove = std::min(2, (int)inClique.size());  // Medium diversification
    } else {
        toRemove = std::min(3, (int)inClique.size());  // Strong diversification late in search
    }

    for (int i = 0; i < toRemove; i++) {
        S.inClique[inClique[i]] = false;
    }
    
    int c = 0;
    for (bool b : S.inClique) {
        if (b) c++;
    }
    S.size = c;
}

Solution TabuSearchMotif::run() {
    initialize();
    Solution S = bestSol;

    if (bestFoundSize == numSequences) {
        return bestSol;
    }

    int iterationCount = 0;
    int stagnationCount = 0;
    Solution restartSolution = S;
    int bestScoreSeen = evaluateMotifSolution(S);

    while (itersSinceImprovement < MaxIter) {
        iterationCount++;
        diversificationCounter++;

        // Check for search stagnation
        if (diversificationCounter >= diversificationInterval) {
            diversify(S);
            diversificationCounter = 0;
            stagnationCount++;
            
            // If we're really stuck, try restarting from the best known solution
            if (stagnationCount >= 5) {
                S = bestSol;
                stagnationCount = 0;
            }
        }

        std::vector<int> C_S = computeC(S);
        std::vector<Solution> N_plus;
        N_plus.reserve(C_S.size());
        
        // Generate augmenting moves
        for (int u : C_S) {
            int seqIdx = G.vertexInfo[u].sequenceIndex;
            bool used = false;
            for (int v = 0; v < G.n; v++) {
                if (S.inClique[v] && G.vertexInfo[v].sequenceIndex == seqIdx) {
                    used = true;
                    break;
                }
            }
            if (used) continue;

            Solution Sprime = S;
            Sprime.inClique[u] = true;
            Sprime.size = S.size + 1;
            N_plus.push_back(Sprime);
        }

        // Generate removing moves
        std::vector<Solution> N_minus;
        {
            std::vector<int> sverts = solutionVertices(S);
            N_minus.reserve(sverts.size());
            for (int v : sverts) {
                Solution Sprime = S;
                Sprime.inClique[v] = false;
                Sprime.size = S.size - 1;
                N_minus.push_back(Sprime);
            }
        }

        // Try augmenting moves first
        Solution chosen = selectBestNeighbor(S, N_plus);
        if (chosen.size == -1) {
            chosen = selectBestNeighbor(S, N_minus);
        }
        
        if (chosen.size == -1) {
            // No valid moves found, diversify
            diversify(S);
            continue;
        }

        bool augmenting = (chosen.size > S.size);
        int changedV = -1;
        
        if (augmenting) {
            for (int i = 0; i < G.n; i++) {
                if (chosen.inClique[i] && !S.inClique[i]) {
                    changedV = i;
                    break;
                }
            }
        } else {
            for (int i = 0; i < G.n; i++) {
                if (S.inClique[i] && !chosen.inClique[i]) {
                    changedV = i;
                    break;
                }
            }
        }

        S = chosen;
        
        // Check if this solution improves upon the best score seen
        int currentScore = evaluateMotifSolution(S);
        if (currentScore > bestScoreSeen) {
            bestScoreSeen = currentScore;
            restartSolution = S;
        }

        updateBestIfNeeded(S);
        updateTabuListAfterMove(S, augmenting, changedV);

        if (bestFoundSize == numSequences) {
            break;
        }
    }

    return bestSol;
}
