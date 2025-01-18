#include "../include/TabuSearchDeterministic.h"
#include <algorithm>
#include <vector>

// Heurystyka Johnsona
static Solution johnsonGreedyHeuristic(const Graph &G) {
    Solution S;
    S.inClique.resize(G.n,false);
    S.size = 0;

    while(true) {
        // compute candidate set C(S)
        std::vector<int> existing;
        for (int v = 0; v < G.n; v++) {
            if (S.inClique[v]) existing.push_back(v);
        }

        std::vector<bool> can(G.n, true);
        for (int v : existing) can[v] = false;

        std::vector<int> C_S;
        for (int v = 0; v < G.n; v++) {
            if (!S.inClique[v]) {
                bool allAdj = true;
                for (int w : existing) {
                    if (!G.adj[w][v]) {
                        allAdj=false; 
                        break;
                    }
                }
                if (allAdj) C_S.push_back(v);
            }
        }
        if (C_S.empty()) break;

        // pick vertex from C(S) that yields largest C(S')
        int bestVal = -1;
        int chosen = -1;
        for (int u : C_S) {
            // hypothetical addition
            std::vector<bool> temp = S.inClique;
            temp[u] = true;
            std::vector<int> newClique;
            for (int i = 0; i < G.n; i++) {
                if (temp[i]) newClique.push_back(i);
            }
            // compute C(S')
            std::vector<bool> can2(G.n, true);
            for (int vv : newClique) can2[vv] = false;
            for (int vv : newClique) {
                for (int x = 0; x < G.n; x++) {
                    if (can2[x] && !G.adj[vv][x]) {
                        can2[x] = false;
                    }
                }
            }
            int countC = 0;
            for (int x = 0; x < G.n; x++) {
                if (can2[x]) countC++;
            }
            if (countC > bestVal) {
                bestVal = countC;
                chosen = u;
            }
        }
        if (chosen == -1) break;
        S.inClique[chosen]=true;
        S.size++;
    }
    return S;
}

TabuSearchDeterministic::TabuSearchDeterministic(const Graph &G, int T1_sz, int T2_sz, int maxIter)
    : TabuSearchBase(G, T1_sz, T2_sz, maxIter) {}

void TabuSearchDeterministic::initialize() {
    bestSol = johnsonGreedyHeuristic(G);
    bestFoundSize = bestSol.size;
    itersSinceImprovement = 0;

    T1_list.clear();
    T1_set.clear();
    T2_list.clear();
    T2_set.clear();

    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

Solution TabuSearchDeterministic::selectBestNeighbor(const Solution &S,
                                                     const std::vector<Solution> &neighbors) {
    int bestVal = -1;
    Solution bestNeighbor;
    bestNeighbor.size = -1;

    for (auto &Sprime : neighbors) {
        // T1
        if (T1_set.find(Sprime) != T1_set.end()) {
            continue;
        }
        bool augmenting = (Sprime.size > S.size);
        if (augmenting) {
            // znajdź wierzchołek dodany
            int addedV = -1;
            for (int i = 0; i < (int)Sprime.inClique.size(); i++) {
                if (Sprime.inClique[i] && !S.inClique[i]) {
                    addedV = i;
                    break;
                }
            }
            if (addedV != -1 && T2_set.find(addedV) != T2_set.end()) {
                // jeśli nie poprawia bestFoundSize to tabu
                if (Sprime.size <= bestFoundSize) {
                    continue;
                }
            }
        }
        // tie-break: liczymy C(S')
        int val = (int)computeC(Sprime).size();
        if (val > bestVal) {
            bestVal = val;
            bestNeighbor = Sprime;
        }
    }
    return bestNeighbor;
}

Solution TabuSearchDeterministic::run() {
    initialize();
    Solution S = bestSol;

    while (itersSinceImprovement < MaxIter) {
        std::vector<int> C_S = computeC(S);

        // N+(S)
        std::vector<Solution> N_plus;
        for (int u : C_S) {
            Solution Sprime = S;
            Sprime.inClique[u] = true;
            Sprime.size = S.size + 1;
            N_plus.push_back(Sprime);
        }

        // N-(S)
        std::vector<Solution> N_minus;
        {
            std::vector<int> sverts = solutionVertices(S);
            for (int v : sverts) {
                Solution Sprime = S;
                Sprime.inClique[v] = false;
                Sprime.size = S.size - 1;
                N_minus.push_back(Sprime);
            }
        }

        // prefer augmenting
        Solution chosen = selectBestNeighbor(S, N_plus);
        if (chosen.size == -1) {
            chosen = selectBestNeighbor(S, N_minus);
        }
        if (chosen.size == -1) {
            // brak ruchu
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
    }
    return bestSol;
}
