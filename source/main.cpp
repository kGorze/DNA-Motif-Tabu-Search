/*
 * Created by konrad_guest on 03/12/2024.
 */


// This code implements a deterministic tabu search approach for the Maximum Clique Problem (MCP)
// as described in Section 3 of the article excerpt provided.
//
// This snippet focuses solely on the deterministic tabu search framework described, 
// i.e., no probabilistic elements and no sampling. It uses the neighborhood structure, 
// tabu lists, aspiration mechanism, and tie-breaking rules as described.
//
// For clarity, we assume the existence of a simple, undirected graph G = (V,E) given by
// an adjacency matrix or adjacency lists. We focus on the essential TS components.
// The code is illustrative and may require further integration and optimization for production use.
//
// The key concepts implemented here include:
// - A solution S is a set of vertices representing a complete subgraph (clique).
// - Neighborhood structure N(S) = N-(S) ∪ N+(S):
//   * N+(S): solutions obtained by adding one vertex u ∈ C(S), where C(S) is the set of vertices adjacent to all in S.
//   * N-(S): solutions obtained by removing one vertex from S.
// - The primary move is always to try to add a vertex (augmenting move) unless it's impossible or tabu-restricted.
// - Tie-breaking rule: choose the neighbor with the maximum |C(S')| (the largest candidate set of further augmentations).
// - Two tabu lists: 
//   * T1: last |T1| visited solutions.
//   * T2: last |T2| vertices removed, but T2 only restricts augmenting moves.
// - Aspiration condition: If a move yields a larger clique than ever found (improves best), it overrides T2 restrictions.
// - If stuck with no improving augmentations, attempt a redirect by considering deletions.
// - Stop after MaxIter steps without improvement.
//
// Inputs:
// - Graph adjacency information
// - Parameters: |T1|, |T2|, MaxIter
//
// Outputs:
// - Best clique found.
//
// Note: This is a core framework and does not include all the auxiliary functions like reading a graph,
// initializing data structures, or advanced memory management. Some simplifications were made for brevity.

#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <deque>

// This code implements a deterministic tabu search approach for the Maximum Clique Problem (MCP)
// as described in Section 3 of the article excerpt provided.
//
// This snippet focuses solely on the deterministic tabu search framework described, 
// i.e., no probabilistic elements and no sampling. It uses the neighborhood structure, 
// tabu lists, aspiration mechanism, and tie-breaking rules as described.
//
// For clarity, we assume the existence of a simple, undirected graph G = (V,E) given by
// an adjacency matrix or adjacency lists. We focus on the essential TS components.
// The code is illustrative and may require further integration and optimization for production use.
//
// The key concepts implemented here include:
// - A solution S is a set of vertices representing a complete subgraph (clique).
// - Neighborhood structure N(S) = N-(S) ∪ N+(S):
//   * N+(S): solutions obtained by adding one vertex u ∈ C(S), where C(S) is the set of vertices adjacent to all in S.
//   * N-(S): solutions obtained by removing one vertex from S.
// - The primary move is always to try to add a vertex (augmenting move) unless it's impossible or tabu-restricted.
// - Tie-breaking rule: choose the neighbor with the maximum |C(S')| (the largest candidate set of further augmentations).
// - Two tabu lists: 
//   * T1: last |T1| visited solutions.
//   * T2: last |T2| vertices removed, but T2 only restricts augmenting moves.
// - Aspiration condition: If a move yields a larger clique than ever found (improves best), it overrides T2 restrictions.
// - If stuck with no improving augmentations, attempt a redirect by considering deletions.
// - Stop after MaxIter steps without improvement.
//
// Inputs:
// - Graph adjacency information
// - Parameters: |T1|, |T2|, MaxIter
//
// Outputs:
// - Best clique found.
//
// Note: This is a core framework and does not include all the auxiliary functions like reading a graph,
// initializing data structures, or advanced memory management. Some simplifications were made for brevity.

#include <iostream>
#include "../include/Graph.h"
#include "../include/TabuSearchDeterministic.h"
#include "../include/TabuSearchProbabilistic.h"
#include "../include/Stabulus.h"
#include "../include/IteratedStabulus.h"
#include "../include/GraphGenerator.h"

int main() {
    // Example usage
    int n = 300;
    double a = 0.95;
    double b = 0.35;

    // Generate a graph with n=100, a=0.25, b=0.75
    Graph G = GraphGenerator::generate(n,a,b);

    // Just print some information about the generated graph
    int edgeCount = 0;
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) {
            if (G.adj[i][j]) edgeCount++;
        }
    }

    std::cout << "Generated graph with " << n << " vertices and " << edgeCount << " edges.\n";

    // Deterministic TS
    {
        TabuSearchDeterministic tsd(G,10,0,1000);
        Solution best=tsd.run();
        std::cout<<"Deterministic TS best clique size: "<<best.size<<"\n";
    }

    // Probabilistic TS
    {
        TabuSearchProbabilistic tsp(G,10,5,1000,2,3);
        Solution best=tsp.run();
        std::cout<<"Probabilistic TS best clique size: "<<best.size<<"\n";
    }

    // Stabulus: find stable set of size k in G
    {
        int k=3; // guess a stable set size
        Stabulus stab(G,k);
        std::vector<int> stableSet;
        bool found=stab.findStableSetOfSizeK(stableSet);
        if(found) {
            std::cout<<"Found stable set of size "<<k<<": ";
            for(int v:stableSet) std::cout<<v<<" ";
            std::cout<<"\n";
        } else {
            std::cout<<"No stable set of size "<<k<<" found.\n";
        }
    }

    // IteratedStabulus: tries to find max stable set by iterating k
    {
        IteratedStabulus istab(G,1,5); // start at k=1, allow 5 fails
        std::vector<int> bestSet;
        int maxSize=istab.runFindMaxStableSet(bestSet);
        std::cout<<"Iterated-Stabulus max stable set size: "<<maxSize<<"\n";
        std::cout<<"Stable set: ";
        for(int v:bestSet) std::cout<<v<<" ";
        std::cout<<"\n";
    }

    return 0;
}
