//
// Created by konrad_guest on 16/12/2024.
//

#include "../include/TabuSearchProbabilistic.h"
#include <random>
#include <chrono>
#include <algorithm>

static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// Reuse Johnson heuristic from deterministic as a starting solution
static Solution johnsonGreedyHeuristicProb(const Graph &G) {
    // same as in deterministic
    Solution S; 
    S.inClique.resize(G.n, false);
    S.size=0;

    while (true) {
        // identical code as above for Johnson
        std::vector<int> C_S;
        {
            std::vector<bool> can(G.n,true);
            for (int v=0; v<G.n; v++) if (S.inClique[v]) can[v]=false;
            std::vector<int> Sverts;
            for (int v=0; v<G.n; v++) if (S.inClique[v]) Sverts.push_back(v);

            for (int v=0; v<G.n; v++) {
                if(!S.inClique[v]) {
                    bool allAdj=true;
                    for (int w : Sverts) {
                        if(!G.adj[w][v]) {allAdj=false;break;}
                    }
                    if(allAdj) C_S.push_back(v);
                }
            }
        }
        if(C_S.empty()) break;

        int bestVal=-1, chosen=-1;
        for (int u : C_S) {
            Solution Sprime=S; Sprime.inClique[u]=true; Sprime.size=S.size+1;

            std::vector<int> C_Sprime;
            {
                std::vector<int> Sverts;
                for (int v=0; v<G.n; v++) if(Sprime.inClique[v]) Sverts.push_back(v);
                std::vector<bool> can(G.n,true);
                for (int v=0; v<G.n; v++) if(Sprime.inClique[v]) can[v]=false;
                for (int w : Sverts) {
                    for (int v=0; v<G.n; v++) {
                        if(can[v]&&!G.adj[w][v]) can[v]=false;
                    }
                }
                for(int v=0; v<G.n; v++) {
                    if(can[v]&&!Sprime.inClique[v]) C_Sprime.push_back(v);
                }
            }
            int val=(int)C_Sprime.size();
            if(val>bestVal) {bestVal=val;chosen=u;}
        }
        if(chosen==-1) break;
        S.inClique[chosen]=true; S.size++;
    }
    return S;
}


TabuSearchProbabilistic::TabuSearchProbabilistic(const Graph &G, int T1_sz, int T2_sz, int maxIter, int k_, int L_)
: TabuSearchBase(G,T1_sz,T2_sz,maxIter), k(k_), L(L_) {}

void TabuSearchProbabilistic::initialize() {
    bestSol = johnsonGreedyHeuristicProb(G);
    bestFoundSize = bestSol.size;
    itersSinceImprovement = 0;
    T1_list.clear(); T1_set.clear(); 
    T2_list.clear(); T2_set.clear();
    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

void TabuSearchProbabilistic::randomShakeUp(Solution &S, int kk) {
    std::vector<int> Sverts = solutionVertices(S);
    if(Sverts.empty()) return;
    int toRemove=std::min(kk,(int)Sverts.size());
    std::shuffle(Sverts.begin(), Sverts.end(), rng);
    for (int i=0; i<toRemove; i++) {
        S.inClique[Sverts[i]]=false;
    }
    int count=0;
    for (bool b : S.inClique) if(b) count++;
    S.size=count;
}

std::vector<Solution> TabuSearchProbabilistic::sampleAugmentingNeighbors(const std::vector<Solution> &N_plus, int LL) {
    std::vector<Solution> candidates=N_plus;
    if ((int)candidates.size()>LL) {
        std::shuffle(candidates.begin(), candidates.end(), rng);
        candidates.resize(LL);
    }
    return candidates;
}

Solution TabuSearchProbabilistic::selectBestNeighbor(const Solution &S,
                                                     const std::vector<Solution> &neighbors) {
    // This is a generic selection without sampling; actual sampling is done in run().
    int bestVal=-1;
    Solution bestNeighbor; bestNeighbor.size=-1;

    for (auto &Sprime : neighbors) {
        if(T1_set.find(Sprime)!=T1_set.end()) continue;
        bool augmenting=(Sprime.size>S.size);
        int addedVertex=-1;
        if(augmenting) {
            for (int i=0; i<(int)Sprime.inClique.size(); i++) {
                if (Sprime.inClique[i]&&!S.inClique[i]) {addedVertex=i;break;}
            }
            if (addedVertex!=-1 && T2_set.find(addedVertex)!=T2_set.end()) {
                if (Sprime.size<=bestFoundSize) continue; 
            }
        }
        int val=(int)computeC(Sprime).size();
        if(val>bestVal) {bestVal=val;bestNeighbor=Sprime;}
    }

    return bestNeighbor;
}


Solution TabuSearchProbabilistic::run() {
    initialize();
    Solution S=bestSol;

    while (itersSinceImprovement<MaxIter) {
        std::vector<int> C_S=computeC(S);
        std::vector<Solution> N_plus;
        for (int u : C_S) {
            Solution Sprime=S; Sprime.inClique[u]=true; Sprime.size=S.size+1;
            N_plus.push_back(Sprime);
        }

        std::vector<Solution> N_minus;
        {
            std::vector<int> Sverts=solutionVertices(S);
            for (int v : Sverts) {
                Solution Sprime=S; Sprime.inClique[v]=false; Sprime.size=S.size-1;
                N_minus.push_back(Sprime);
            }
        }

        if(N_plus.empty()) {
            // local max
            randomShakeUp(S,k);
            updateBestIfNeeded(S);
            updateTabuListAfterMove(S,false,-1);
            continue;
        } else {
            std::vector<Solution> sampleSet=sampleAugmentingNeighbors(N_plus,L);

            int bestVal=-1;
            Solution bestSample; bestSample.size=-1;
            for (auto &Sprime : sampleSet) {
                if(T1_set.find(Sprime)!=T1_set.end()) continue;
                bool augmenting=(Sprime.size>S.size);
                int addedVertex=-1;
                if(augmenting) {
                    for (int i=0; i<(int)Sprime.inClique.size(); i++) {
                        if (Sprime.inClique[i]&&!S.inClique[i]) {addedVertex=i;break;}
                    }
                    if(addedVertex!=-1 && T2_set.find(addedVertex)!=T2_set.end()) {
                        if (Sprime.size<=bestFoundSize) continue;
                    }
                }
                int val=(int)computeC(Sprime).size();
                if(val>bestVal) {bestVal=val;bestSample=Sprime;}
            }

            if(bestSample.size==-1) {
                // no feasible from sample
                randomShakeUp(S,k);
                updateBestIfNeeded(S);
                updateTabuListAfterMove(S,false,-1);
                continue;
            } else {
                // perform move to bestSample
                bool augmenting=(bestSample.size>S.size);
                int changedVertex=-1;
                if(augmenting) {
                    for(int i=0; i<G.n; i++) {
                        if(bestSample.inClique[i]&&!S.inClique[i]) {changedVertex=i;break;}
                    }
                } else {
                    for(int i=0; i<G.n; i++) {
                        if(S.inClique[i]&&!bestSample.inClique[i]) {changedVertex=i;break;}
                    }
                }
                S=bestSample;
                updateBestIfNeeded(S);
                updateTabuListAfterMove(S,augmenting,changedVertex);
            }
        }
    }
    return bestSol;
}

