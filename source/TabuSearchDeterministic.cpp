//
// Created by konrad_guest on 16/12/2024.
//
#include "../include/TabuSearchDeterministic.h"

// Johnson's greedy heuristic
static Solution johnsonGreedyHeuristic(const Graph &G) {
    Solution S; 
    S.inClique.resize(G.n, false);
    S.size = 0;

    while (true) {
        // choose vertex in C(S) that yields largest C(S')
        // from previous code
        std::vector<int> C_S;
        {
            std::vector<bool> can(G.n,true);
            for (int v = 0; v < G.n; v++) if (S.inClique[v]) can[v]=false;
            std::vector<int> Sverts;
            for (int v=0; v<G.n; v++) if (S.inClique[v]) Sverts.push_back(v);

            for (int v=0; v<G.n; v++) {
                if (!S.inClique[v]) {
                    bool allAdj = true;
                    for (int w : Sverts) {
                        if (!G.adj[w][v]) {allAdj=false;break;}
                    }
                    if (allAdj) C_S.push_back(v);
                }
            }
        }
        if (C_S.empty()) break;

        int bestVal = -1, chosen = -1;
        for (int u : C_S) {
            Solution Sprime = S;
            Sprime.inClique[u] = true;
            Sprime.size = S.size+1;

            // compute C(S')
            std::vector<int> C_Sprime;
            {
                std::vector<int> Sverts;
                for (int v=0; v<G.n; v++) if (Sprime.inClique[v]) Sverts.push_back(v);
                std::vector<bool> can(G.n,true);
                for (int v=0; v<G.n; v++) if (Sprime.inClique[v]) can[v]=false;
                for (int w : Sverts) {
                    for (int v=0; v<G.n; v++) {
                        if (can[v] && !G.adj[w][v]) {
                            can[v]=false;
                        }
                    }
                }
                for (int v=0; v<G.n; v++) {
                    if (can[v] && !Sprime.inClique[v]) C_Sprime.push_back(v);
                }
            }

            int val = (int)C_Sprime.size();
            if (val>bestVal) {
                bestVal = val;
                chosen = u;
            }
        }
        if (chosen==-1) break;
        S.inClique[chosen]=true; S.size++;
    }

    return S;
}

TabuSearchDeterministic::TabuSearchDeterministic(const Graph &G, int T1_sz, int T2_sz, int maxIter)
: TabuSearchBase(G,T1_sz,T2_sz,maxIter) {}

void TabuSearchDeterministic::initialize() {
    // start from empty set and perform Johnson's heuristic
    bestSol = johnsonGreedyHeuristic(G);
    bestFoundSize = bestSol.size;
    itersSinceImprovement = 0;

    // Clear tabu lists
    T1_list.clear(); T1_set.clear(); 
    T2_list.clear(); T2_set.clear();

    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

Solution TabuSearchDeterministic::selectBestNeighbor(const Solution &S,
                                                     const std::vector<Solution> &neighbors) {
    int bestVal = -1;
    Solution bestNeighbor; bestNeighbor.size=-1;

    for (auto &Sprime : neighbors) {
        // Check tabu T1
        if (T1_set.find(Sprime)!=T1_set.end()) continue;

        bool augmenting = (Sprime.size>S.size);
        if (augmenting) {
            // Identify added vertex
            int addedVertex=-1;
            for (int i=0; i<(int)Sprime.inClique.size(); i++) {
                if (Sprime.inClique[i]&&!S.inClique[i]) {addedVertex=i; break;}
            }
            if (addedVertex!=-1 && T2_set.find(addedVertex)!=T2_set.end()) {
                // Tabu unless improves best
                if (Sprime.size<=bestFoundSize) continue;
            }
        }

        int val = (int)computeC(Sprime).size();
        if (val>bestVal) {
            bestVal=val;
            bestNeighbor=Sprime;
        }
    }

    return bestNeighbor;
}

Solution TabuSearchDeterministic::run() {
    initialize();
    Solution S = bestSol;

    while (itersSinceImprovement<MaxIter) {
        std::vector<int> C_S = computeC(S);

        // N+(S)
        std::vector<Solution> N_plus;
        for (int u : C_S) {
            Solution Sprime=S; Sprime.inClique[u]=true; Sprime.size=S.size+1;
            N_plus.push_back(Sprime);
        }

        // N-(S)
        std::vector<Solution> N_minus;
        {
            std::vector<int> Sverts=solutionVertices(S);
            for (int v : Sverts) {
                Solution Sprime=S; Sprime.inClique[v]=false; Sprime.size=S.size-1;
                N_minus.push_back(Sprime);
            }
        }

        Solution chosenNeighbor; chosenNeighbor.size=-1;
        if (!N_plus.empty()) {
            chosenNeighbor=selectBestNeighbor(S,N_plus);
        }
        if (chosenNeighbor.size==-1) {
            chosenNeighbor=selectBestNeighbor(S,N_minus);
        }

        if (chosenNeighbor.size==-1) {
            // no moves found
            // Attempt redirection: 
            if (!N_minus.empty()) {
                // arbitrarily choose any from N_minus if all else fails
                chosenNeighbor=N_minus.front();
            } else {
                break; // stuck
            }
        }

        bool augmenting=(chosenNeighbor.size>S.size);
        int changedVertex=-1;
        if (augmenting) {
            for (int i=0; i<G.n; i++) {
                if (chosenNeighbor.inClique[i]&&!S.inClique[i]) {changedVertex=i;break;}
            }
        } else {
            for (int i=0; i<G.n; i++) {
                if (S.inClique[i]&&!chosenNeighbor.inClique[i]) {changedVertex=i;break;}
            }
        }

        S=chosenNeighbor;
        updateBestIfNeeded(S);
        updateTabuListAfterMove(S,augmenting,changedVertex);
    }

    return bestSol;
}

