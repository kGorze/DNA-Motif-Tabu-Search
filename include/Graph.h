//
// Created by konrad_guest on 16/12/2024.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

// This struct holds all the extra metadata needed for motif finding.
struct VertexData {
    int sequenceIndex;   // which DNA sequence
    int position;        // position in that sequence
    std::string kmer;    // the k-mer string
};

class Graph {
public:
    int n; // number of vertices

    // Adjacency matrix: adj[u][v] = true if edge (u,v) exists
    std::vector<std::vector<bool>> adj;

    // Additional data for each vertex, used in motif context
    std::vector<VertexData> vertexInfo;

    // Constructor for an empty graph with n vertices
    // (Will fill adjacency with false, and you can fill vertexInfo externally)
    Graph(int n_);

    // Add an undirected edge (u,v)
    void add_edge(int u, int v);

    // Check if there is an edge (u,v)
    bool is_edge(int u, int v) const;

    // For convenience in motif-related tasks: check if two vertices can be connected
    // i.e., identical k-mers, different sequences, position difference constraint, etc.
    // Not used in the old code, but can be used in new code if needed.
    static bool canConnect(const VertexData &a, const VertexData &b, int positionMultiplier);
};
#endif //GRAPH_H
