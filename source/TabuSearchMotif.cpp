//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/TabuSearchMotif.h"
#include <algorithm>

// A simple heuristic start: pick any valid vertex from each sequence if possible
static Solution naiveMotifInitialization(const Graph &G, int numSequences) {
    // We'll pick the first vertex we find for each sequence, if feasible
    // This is a naive approach
    Solution S;
    S.inClique.resize(G.n, false);
    S.size = 0;

    std::vector<bool> sequenceUsed(numSequences, false);

    for (int v = 0; v < G.n; v++) {
        int seqIndex = G.vertexInfo[v].sequenceIndex;
        if (!sequenceUsed[seqIndex]) {
            // check if it is adjacent to all already chosen vertices
            bool canPick = true;
            for (int u = 0; u < G.n; u++) {
                if (S.inClique[u]) {
                    if (!G.adj[u][v]) {
                        canPick = false;
                        break;
                    }
                }
            }
            if (canPick) {
                S.inClique[v] = true;
                S.size++;
                sequenceUsed[seqIndex] = true;
                if (S.size == numSequences) break;
            }
        }
    }

    return S;
}

TabuSearchMotif::TabuSearchMotif(const Graph &G, int T1_sz, int T2_sz, int maxIter, int numSeq)
    : TabuSearchBase(G, T1_sz, T2_sz, maxIter), numSequences(numSeq) {}

void TabuSearchMotif::initialize() {
    bestSol = naiveMotifInitialization(G, numSequences);
    bestFoundSize = bestSol.size;
    itersSinceImprovement = 0;

    T1_list.clear();
    T1_set.clear();
    T2_list.clear();
    T2_set.clear();

    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

bool TabuSearchMotif::validMotifSolution(const Solution &S) const {
    // Check for no repeated sequence
    // Then check it is a clique
    std::vector<int> used(numSequences, 0);

    // also collect vertices
    std::vector<int> verts;
    for (int i = 0; i < G.n; i++) {
        if (S.inClique[i]) {
            int seqIdx = G.vertexInfo[i].sequenceIndex;
            used[seqIdx]++;
            // If any sequence used more than once => invalid
            if (used[seqIdx] > 1) return false;
            verts.push_back(i);
        }
    }

    // clique check
    for (size_t i = 0; i < verts.size(); i++) {
        for (size_t j = i + 1; j < verts.size(); j++) {
            if (!G.adj[verts[i]][verts[j]]) {
                return false;
            }
        }
    }

    return true;
}

Solution TabuSearchMotif::selectBestNeighbor(const Solution &S, const std::vector<Solution> &neighbors) {
    int bestVal = -1;
    Solution bestNeighbor;
    bestNeighbor.size = -1;

    for (auto &Sprime : neighbors) {
        if (!validMotifSolution(Sprime)) {
            continue; // must have exactly 0 or 1 vertex from each sequence
        }
        if (T1_set.find(Sprime) != T1_set.end()) {
            continue; // taboo
        }
        bool augmenting = (Sprime.size > S.size);
        if (augmenting) {
            // find added vertex
            int addedV = -1;
            for (int i = 0; i < G.n; i++) {
                if (Sprime.inClique[i] && !S.inClique[i]) {
                    addedV = i;
                    break;
                }
            }
            if (addedV != -1 && T2_set.find(addedV) != T2_set.end()) {
                // only allow if it improves best
                if (Sprime.size <= bestFoundSize) {
                    continue;
                }
            }
        }
        // tie-break: prefer bigger clique, then bigger C(S')
        // (or just bigger C(S') if size is same)
        if (Sprime.size > bestVal) {
            bestVal = Sprime.size;
            bestNeighbor = Sprime;
        } else if (Sprime.size == bestVal) {
            // secondary tie-break by size of C(Sprime)
            int csize1 = (int)computeC(Sprime).size();
            if (bestNeighbor.size >= 0) {
                int csize2 = (int)computeC(bestNeighbor).size();
                if (csize1 > csize2) {
                    bestNeighbor = Sprime;
                }
            } else {
                bestNeighbor = Sprime;
            }
        }
    }
    return bestNeighbor;
}

Solution TabuSearchMotif::run() {
    initialize();
    Solution S = bestSol;

    // If we already have the perfect clique of size = numSequences, we can stop
    if (bestFoundSize == numSequences) {
        return bestSol;
    }

    while (itersSinceImprovement < MaxIter) {
        // Build neighbors
        std::vector<int> C_S = computeC(S);

        // N+(S)
        std::vector<Solution> N_plus;
        for (int u : C_S) {
            // but only add if that sequence isn't used
            int seqIdx = G.vertexInfo[u].sequenceIndex;
            // check if S already uses that seq
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

        // N-(S)
        std::vector<Solution> N_minus;
        {
            std::vector<int> Sverts = solutionVertices(S);
            for (int v : Sverts) {
                Solution Sprime = S;
                Sprime.inClique[v] = false;
                Sprime.size = S.size - 1;
                N_minus.push_back(Sprime);
            }
        }

        // prefer augmenting
        Solution chosen = selectBestNeighbor(S, N_plus);

        // if none chosen, try removing
        if (chosen.size == -1) {
            chosen = selectBestNeighbor(S, N_minus);
        }

        // if still none found, stuck
        if (chosen.size == -1) {
            // no moves
            break;
        }

        bool augmenting = (chosen.size > S.size);
        int changedVertex = -1;
        if (augmenting) {
            for (int i = 0; i < G.n; i++) {
                if (chosen.inClique[i] && !S.inClique[i]) {
                    changedVertex = i;
                    break;
                }
            }
        } else {
            for (int i = 0; i < G.n; i++) {
                if (S.inClique[i] && !chosen.inClique[i]) {
                    changedVertex = i;
                    break;
                }
            }
        }

        S = chosen;
        updateBestIfNeeded(S);
        updateTabuListAfterMove(S, augmenting, changedVertex);

        if (bestFoundSize == numSequences) {
            // We have a perfect solution
            break;
        }
    }

    return bestSol;
}
