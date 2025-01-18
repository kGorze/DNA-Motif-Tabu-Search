#include "../include/TabuSearchDeterministic.h"
#include <vector>
#include <algorithm>

// Heurystyka Johnsona do wygenerowania kliki startowej
static Solution johnsonGreedyHeuristic(const Graph &G) {
    Solution S;
    S.inClique.resize(G.n, false);
    S.size = 0;

    while (true) {
        // Wyznaczamy C(S)
        std::vector<int> C_S;
        {
            std::vector<int> existing;
            for (int v = 0; v < G.n; v++)
                if (S.inClique[v]) existing.push_back(v);

            std::vector<bool> can(G.n, true);
            for (int v : existing) can[v] = false;

            for (int v = 0; v < G.n; v++) {
                if (!S.inClique[v]) {
                    bool allAdj = true;
                    for (int w : existing) {
                        if (!G.adj[w][v]) { allAdj = false; break; }
                    }
                    if (!allAdj) can[v] = false;
                }
            }
            for (int v = 0; v < G.n; v++) {
                if (can[v] && !S.inClique[v]) {
                    C_S.push_back(v);
                }
            }
        }
        if (C_S.empty()) break;

        // Wybieramy wierzchołek z C(S) maksymalizujący |C(S')|
        int bestVal = -1;
        int chosen = -1;
        for (int u : C_S) {
            std::vector<bool> temp = S.inClique;
            temp[u] = true;
            std::vector<int> newCliqueVerts;
            for (int i = 0; i < G.n; i++) {
                if (temp[i]) newCliqueVerts.push_back(i);
            }
            // obliczamy C(S')
            std::vector<bool> can(G.n, true);
            for (int vv : newCliqueVerts) can[vv] = false;
            for (int vv : newCliqueVerts) {
                for (int x = 0; x < G.n; x++) {
                    if (can[x] && !G.adj[vv][x]) {
                        can[x] = false;
                    }
                }
            }
            int countC = 0;
            for (int x = 0; x < G.n; x++) {
                if (can[x]) countC++;
            }
            if (countC > bestVal) {
                bestVal = countC;
                chosen = u;
            }
        }

        if (chosen == -1) break;
        S.inClique[chosen] = true;
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
        // sprawdzenie tabu T1
        if (T1_set.find(Sprime) != T1_set.end()) {
            continue;
        }
        bool augmenting = (Sprime.size > S.size);

        // jeśli augmenting, sprawdzamy T2
        if (augmenting) {
            int addedV = -1;
            for (int i = 0; i < (int)Sprime.inClique.size(); i++) {
                if (Sprime.inClique[i] && !S.inClique[i]) {
                    addedV = i;
                    break;
                }
            }
            // jeśli wierzchołek w T2 i nie poprawiamy bestFoundSize, pomijamy
            if (addedV != -1 && T2_set.find(addedV) != T2_set.end()) {
                if (Sprime.size <= bestFoundSize) {
                    continue;
                }
            }
        }

        // Kryterium wyboru sąsiada: maksymalizujemy |C(S')|
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
            std::vector<int> Sverts = solutionVertices(S);
            for (int v : Sverts) {
                Solution Sprime = S;
                Sprime.inClique[v] = false;
                Sprime.size = S.size - 1;
                N_minus.push_back(Sprime);
            }
        }

        // Najpierw próbujemy augmentacji
        Solution chosen = selectBestNeighbor(S, N_plus);

        // Jeśli nie ma augmentacji, próbujemy usunięć
        if (chosen.size == -1) {
            chosen = selectBestNeighbor(S, N_minus);
        }

        // Jeśli dalej nic, koniec
        if (chosen.size == -1) {
            if (!N_minus.empty()) {
                chosen = N_minus.front(); 
            } else {
                break;
            }
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
