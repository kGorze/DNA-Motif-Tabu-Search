#include "../include/TabuSearchProbabilistic.h"
#include <random>
#include <chrono>
#include <algorithm>

static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

static Solution johnsonGreedyHeuristicProb(const Graph &G) {
    // identyczne jak w deterministic, do uproszczenia
    Solution S;
    S.inClique.resize(G.n,false);
    S.size=0;

    while(true) {
        std::vector<int> existing;
        for (int v=0; v<G.n; v++){
            if(S.inClique[v]) existing.push_back(v);
        }
        std::vector<bool> can(G.n,true);
        for (int v:existing) can[v]=false;

        std::vector<int> C_S;
        for(int v=0; v<G.n; v++){
            if(!S.inClique[v]){
                bool allAdj=true;
                for (int w:existing){
                    if(!G.adj[w][v]){allAdj=false;break;}
                }
                if(allAdj) C_S.push_back(v);
            }
        }
        if(C_S.empty()) break;

        int bestVal=-1, chosen=-1;
        for(int u: C_S){
            std::vector<bool> temp=S.inClique;
            temp[u]=true;
            std::vector<int> newClique;
            for(int i=0; i<G.n; i++){
                if(temp[i]) newClique.push_back(i);
            }
            std::vector<bool> can2(G.n,true);
            for(int vv: newClique) can2[vv]=false;
            for(int vv: newClique){
                for(int x=0; x<G.n; x++){
                    if(can2[x] && !G.adj[vv][x]){
                        can2[x]=false;
                    }
                }
            }
            int countC=0;
            for(int x=0; x<G.n; x++){
                if(can2[x]) countC++;
            }
            if(countC>bestVal){
                bestVal=countC;
                chosen=u;
            }
        }
        if(chosen==-1) break;
        S.inClique[chosen]=true;
        S.size++;
    }
    return S;
}

TabuSearchProbabilistic::TabuSearchProbabilistic(const Graph &G, int T1_sz, int T2_sz, int maxIter, int k_, int L_)
: TabuSearchBase(G,T1_sz,T2_sz,maxIter), k(k_), L(L_) {}

void TabuSearchProbabilistic::initialize() {
    bestSol = johnsonGreedyHeuristicProb(G);
    bestFoundSize = bestSol.size;
    itersSinceImprovement=0;
    T1_list.clear();
    T1_set.clear();
    T2_list.clear();
    T2_set.clear();
    T1_list.push_back(bestSol);
    T1_set.insert(bestSol);
}

void TabuSearchProbabilistic::randomShakeUp(Solution &S, int kk){
    std::vector<int> sverts;
    for (int i=0; i<(int)S.inClique.size(); i++){
        if(S.inClique[i]) sverts.push_back(i);
    }
    if(sverts.empty()) return;
    std::shuffle(sverts.begin(), sverts.end(), rng);

    int toRemove=std::min(kk,(int)sverts.size());
    for(int i=0; i<toRemove; i++){
        S.inClique[sverts[i]] = false;
    }
    int c=0;
    for(bool b: S.inClique) if(b) c++;
    S.size=c;
}

std::vector<Solution> TabuSearchProbabilistic::sampleAugmentingNeighbors(const std::vector<Solution> &N_plus, int LL){
    std::vector<Solution> sample=N_plus;
    if((int)sample.size() > LL){
        std::shuffle(sample.begin(), sample.end(), rng);
        sample.resize(LL);
    }
    return sample;
}

Solution TabuSearchProbabilistic::selectBestNeighbor(const Solution &S,
                                                     const std::vector<Solution> &neighbors){
    int bestVal=-1;
    Solution bestNeighbor; 
    bestNeighbor.size=-1;

    for(auto &Sprime : neighbors){
        if(T1_set.find(Sprime)!=T1_set.end()) continue;
        bool augmenting=(Sprime.size>S.size);
        if(augmenting){
            int addedV=-1;
            for(int i=0; i<(int)Sprime.inClique.size(); i++){
                if(Sprime.inClique[i] && !S.inClique[i]){addedV=i;break;}
            }
            if(addedV!=-1 && T2_set.find(addedV)!=T2_set.end()){
                if(Sprime.size<=bestFoundSize) continue;
            }
        }
        int val=(int)computeC(Sprime).size();
        if(val>bestVal){
            bestVal=val;
            bestNeighbor=Sprime;
        }
    }
    return bestNeighbor;
}

Solution TabuSearchProbabilistic::run(){
    initialize();
    Solution S=bestSol;

    while(itersSinceImprovement<MaxIter){
        std::vector<int> C_S=computeC(S);

        std::vector<Solution> N_plus;
        for(int u : C_S){
            Solution Sprime=S;
            Sprime.inClique[u]=true;
            Sprime.size=S.size+1;
            N_plus.push_back(Sprime);
        }

        std::vector<Solution> N_minus;
        {
            std::vector<int> sverts;
            for(int i=0; i<(int)S.inClique.size(); i++){
                if(S.inClique[i]) sverts.push_back(i);
            }
            for(int v : sverts){
                Solution Sprime=S;
                Sprime.inClique[v]=false;
                Sprime.size=S.size-1;
                N_minus.push_back(Sprime);
            }
        }

        if(N_plus.empty()){
            randomShakeUp(S,k);
            updateBestIfNeeded(S);
            updateTabuListAfterMove(S,false,-1);
            continue;
        } else {
            std::vector<Solution> sampleSet=sampleAugmentingNeighbors(N_plus,L);
            int bestVal=-1;
            Solution bestSample; 
            bestSample.size=-1;

            for (auto &Sprime : sampleSet){
                if(T1_set.find(Sprime)!=T1_set.end()) continue;
                bool augmenting=(Sprime.size>S.size);
                int addedV=-1;
                if(augmenting){
                    for(int i=0; i<(int)Sprime.inClique.size(); i++){
                        if(Sprime.inClique[i] && !S.inClique[i]){addedV=i;break;}
                    }
                    if(addedV!=-1 && T2_set.find(addedV)!=T2_set.end()){
                        if(Sprime.size<=bestFoundSize) continue;
                    }
                }
                int val=(int)computeC(Sprime).size();
                if(val>bestVal){
                    bestVal=val;
                    bestSample=Sprime;
                }
            }
            if(bestSample.size==-1){
                randomShakeUp(S,k);
                updateBestIfNeeded(S);
                updateTabuListAfterMove(S,false,-1);
                continue;
            } else {
                bool augmenting=(bestSample.size>S.size);
                int changedVertex=-1;
                if(augmenting){
                    for(int i=0; i<G.n; i++){
                        if(bestSample.inClique[i] && !S.inClique[i]){
                            changedVertex=i;break;
                        }
                    }
                } else {
                    for(int i=0; i<G.n; i++){
                        if(S.inClique[i] && !bestSample.inClique[i]){
                            changedVertex=i;break;
                        }
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
