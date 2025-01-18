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

#include <iostream>
#include <vector>
#include <string>
#include "../include/Graph.h"
#include "../include/GraphGenerator.h"
#include "../include/TabuSearchDeterministic.h"
#include "../include/TabuSearchProbabilistic.h"
#include "../include/Stabulus.h"
#include "../include/IteratedStabulus.h"
#include "../include/FastaQualParser.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/TabuSearchMotif.h"

int main(int argc, char** argv) {
    // ------------------------------------------------------------
    // Part 1: Demonstrate existing code on a random graph
    // ------------------------------------------------------------
    std::cout << "=== Maximum Clique via Tabu Search (Random Graph Demo) ===\n";

    int n = 100;
    double a = 0.25;
    double b = 0.75;
    Graph G = GraphGenerator::generate(n, a, b);

    // Count edges just for info
    int edgeCount = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (G.adj[i][j]) edgeCount++;
        }
    }
    std::cout << "Random graph with " << n << " vertices, " << edgeCount << " edges.\n";

    // Deterministic TS for Max Clique
    {
        TabuSearchDeterministic tsd(G, 10, 5, 500);
        Solution best = tsd.run();
        std::cout << "Deterministic TS best clique size: " << best.size << "\n";
    }

    // Probabilistic TS for Max Clique
    {
        TabuSearchProbabilistic tsp(G, 10, 5, 500, 2, 3);
        Solution best = tsp.run();
        std::cout << "Probabilistic TS best clique size: " << best.size << "\n";
    }

    // ------------------------------------------------------------
    // Part 2: Demonstrate motif-finding extension
    // ------------------------------------------------------------
    std::cout << "\n=== DNA Motif Finding Demo ===\n";
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <fasta_file> <qual_file>\n";
        std::cout << "Skipping motif demo due to missing arguments.\n";
        return 0;
    }

    // Example parameters
    int kMin = 5;
    int kMax = 6;
    int qualityThreshold = 20;
    int positionMultiplier = 10;  // for difference <= 10*k

    try {
        // 1) Parse FASTA/QUAL
        std::string fastaFile = argv[1];
        std::string qualFile = argv[2];
        std::vector<DNASequence> sequences = FastaQualParser::parseFastaAndQual(fastaFile, qualFile);
        int seqCount = (int)sequences.size();

        std::cout << "Loaded " << seqCount << " sequences.\n";
        for (int i = 0; i < seqCount; i++) {
            std::cout << "  Seq " << i << ": " << sequences[i].name 
                      << " (length=" << sequences[i].length() << ")\n";
        }

        // 2) Build motif graph
        MotifGraphBuilder builder(sequences, kMin, kMax, qualityThreshold, positionMultiplier);
        Graph motifGraph = builder.build();
        std::cout << "Motif graph built with " << motifGraph.n << " vertices.\n";

        // 3) Run specialized TabuSearch for motif
        // Each input sequence must appear exactly once => we want a clique of size = seqCount
        TabuSearchMotif motifSearch(motifGraph, 10, 5, 500, seqCount);
        Solution motifSol = motifSearch.run();

        // 4) Output the found motif
        if (motifSol.size == seqCount) {
            std::cout << "Found a motif clique of size " << motifSol.size << " (one vertex per sequence).\n";
        } else {
            std::cout << "Found a near-solution of size " << motifSol.size
                      << ", expected " << seqCount << " if perfect.\n";
        }

        // Print the details of each vertex in the found set
        std::vector<int> chosenVerts;
        for (int i = 0; i < motifGraph.n; i++) {
            if (motifSol.inClique[i]) {
                chosenVerts.push_back(i);
            }
        }
        for (int idx : chosenVerts) {
            const VertexData &vd = motifGraph.vertexInfo[idx];
            std::cout << "  Sequence=" << vd.sequenceIndex
                      << " Position=" << vd.position
                      << " K-mer=" << vd.kmer << "\n";
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "ERROR in motif-finding demo: " << ex.what() << "\n";
    }

    return 0;
}

