#ifndef TABUSEARCHBASE_H
#define TABUSEARCHBASE_H

#include <vector>
#include <deque>
#include <unordered_set>
#include "Graph.h"

struct Solution {
    std::vector<bool> inClique;
    int size;
};

struct SolutionHash {
    std::size_t operator()(const Solution &S) const;
};

struct SolutionEq {
    bool operator()(const Solution &A, const Solution &B) const;
};

class TabuSearchBase {
protected:
    const Graph &G;
    int T1_size;
    int T2_size;
    int MaxIter;

    std::deque<Solution> T1_list;
    std::unordered_set<Solution, SolutionHash, SolutionEq> T1_set;
    std::deque<int> T2_list;
    std::unordered_set<int> T2_set;

    int bestFoundSize;
    Solution bestSol;
    int itersSinceImprovement;

    std::vector<int> solutionVertices(const Solution &S) const;
    std::vector<int> computeC(const Solution &S) const;
    Solution makeSolution(const std::vector<int> &verts) const;
    int upperBoundFromS(const Solution &S) const;
    bool isClique(const std::vector<int> &solVec) const;

    virtual void initialize() = 0;
    virtual Solution selectBestNeighbor(const Solution &S,
                                        const std::vector<Solution> &neighbors) = 0;

    void updateBestIfNeeded(const Solution &S);
    void updateTabuListAfterMove(const Solution &S, bool augmenting, int changedVertex);

public:
    TabuSearchBase(const Graph &graph, int T1_sz, int T2_sz, int maxIter);
    virtual ~TabuSearchBase() = default;

    virtual Solution run() = 0;
};

#endif // TABUSEARCHBASE_H
