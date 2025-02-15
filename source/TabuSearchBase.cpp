#include "../include/TabuSearchBase.h"
#include <algorithm>

std::size_t SolutionHash::operator()(const Solution &S) const {
    // Prosty hash polinomialny
    std::size_t h = 0;
    for (bool b : S.inClique) {
        h = h * 131 + (b ? 1 : 0);
    }
    return h;
}

bool SolutionEq::operator()(const Solution &A, const Solution &B) const {
    return A.inClique == B.inClique;
}

TabuSearchBase::TabuSearchBase(const Graph &graph, int T1_sz, int T2_sz, int maxIter)
    : G(graph), T1_size(T1_sz), T2_size(T2_sz), MaxIter(maxIter)
{
    bestFoundSize = 0;
    itersSinceImprovement = 0;
}

std::vector<int> TabuSearchBase::solutionVertices(const Solution &S) const {
    std::vector<int> verts;
    for (int i = 0; i < (int)S.inClique.size(); i++) {
        if (S.inClique[i]) verts.push_back(i);
    }
    return verts;
}

std::vector<int> TabuSearchBase::computeC(const Solution &S) const {
    // C(S): wierzchołki spoza S, które są sąsiadami wszystkich wierzchołków S
    if (S.size == 0) {
        std::vector<int> allV(G.n);
        for (int i = 0; i < G.n; i++) allV[i] = i;
        return allV;
    }

    std::vector<int> inS = solutionVertices(S);
    std::vector<bool> can(G.n, true);

    for (int v : inS) {
        can[v] = false; 
    }

    for (int v : inS) {
        for (int x = 0; x < G.n; x++) {
            if (can[x] && !G.is_edge(v,x)) {
                can[x] = false;
            }
        }
    }

    std::vector<int> result;
    for (int x = 0; x < G.n; x++) {
        if (can[x]) result.push_back(x);
    }
    return result;
}

Solution TabuSearchBase::makeSolution(const std::vector<int> &verts) const {
    Solution S;
    S.inClique.resize(G.n, false);
    for (int v : verts) {
        S.inClique[v] = true;
    }
    S.size = (int)verts.size();
    return S;
}

int TabuSearchBase::upperBoundFromS(const Solution &S) const {
    return S.size + (int)computeC(S).size();
}

bool TabuSearchBase::isClique(const std::vector<int> &solVec) const {
    for (size_t i = 0; i < solVec.size(); i++) {
        for (size_t j = i+1; j < solVec.size(); j++) {
            if (!G.is_edge(solVec[i], solVec[j])) return false;
        }
    }
    return true;
}

void TabuSearchBase::updateBestIfNeeded(const Solution &S) {
    if (S.size > bestFoundSize) {
        bestFoundSize = S.size;
        bestSol = S;
        itersSinceImprovement = 0;
    } else {
        itersSinceImprovement++;
    }
}

void TabuSearchBase::updateTabuListAfterMove(const Solution &S, bool augmenting, int changedVertex) {
    // T1
    T1_list.push_back(S);
    T1_set.insert(S);
    if ((int)T1_list.size() > T1_size) {
        T1_set.erase(T1_list.front());
        T1_list.pop_front();
    }

    // T2 - jeśli ruch to usunięcie
    if (!augmenting && changedVertex != -1) {
        T2_list.push_back(changedVertex);
        T2_set.insert(changedVertex);
        if ((int)T2_list.size() > T2_size) {
            int oldv = T2_list.front();
            T2_list.pop_front();
            T2_set.erase(oldv);
        }
    }
}
