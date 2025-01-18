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

    std::vector<bool> used(numSeq, false);

    for (int v = 0; v < G.n; v++) {
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
      diversificationInterval(50),
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
    for (size_t i = 0; i < verts.size(); i++) {
        for (size_t j = i+1; j < verts.size(); j++) {
            if (!G.is_edge(verts[i], verts[j])) {
                return false;
            }
        }
    }
    return true;
}

int TabuSearchMotif::evaluateMotifSolution(const Solution& S) const {
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
        for (size_t j = i+1; j < vertices.size(); j++) {
            const auto &k1 = G.vertexInfo[vertices[i]].kmer;
            const auto &k2 = G.vertexInfo[vertices[j]].kmer;
            totalDiff += countDifferences(k1, k2);
            pairs++;
        }
    }
    if (pairs > 0) {
        double avgDiff = (double)totalDiff / pairs;
        score -= (int)(avgDiff * 10.0); 
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
    int toRemove = std::min((int)inClique.size(), 2);
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

    while (itersSinceImprovement < MaxIter) {
        iterationCount++;
        diversificationCounter++;

        if (diversificationCounter >= diversificationInterval) {
            diversify(S);
            diversificationCounter = 0;
        }

        std::vector<int> C_S = computeC(S);
        std::vector<Solution> N_plus;
        N_plus.reserve(C_S.size());
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

        Solution chosen = selectBestNeighbor(S, N_plus);
        if (chosen.size == -1) {
            chosen = selectBestNeighbor(S, N_minus);
        }
        if (chosen.size == -1) {
            break;
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
        updateBestIfNeeded(S);
        updateTabuListAfterMove(S, augmenting, changedV);

        if (bestFoundSize == numSequences) {
            break;
        }
    }

    return bestSol;
}
