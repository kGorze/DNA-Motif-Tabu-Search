#include "../include/Graph.h"

Graph::Graph(int n_) : n(n_), adj(n_, std::vector<bool>(n_, false)), vertexInfo(n_) {}

void Graph::add_edge(int u, int v) {
    // Avoid self-loops or out-of-range
    if (u == v) return;
    if (u < 0 || v < 0 || u >= n || v >= n) return;

    adj[u][v] = true;
    adj[v][u] = true;
}

bool Graph::is_edge(int u, int v) const {
    if (u < 0 || v < 0 || u >= n || v >= n) return false;
    return adj[u][v];
}

bool Graph::canConnect(const VertexData &a, const VertexData &b, int positionMultiplier) {
    // Example rule for edges in motif context:
    // 1) k-mers must be identical
    // 2) from different sequences
    // 3) position difference <= positionMultiplier * k
    if (a.sequenceIndex == b.sequenceIndex) return false;
    if (a.kmer != b.kmer) return false;
    int k = (int)a.kmer.size();
    int diff = (a.position > b.position)
               ? (a.position - b.position)
               : (b.position - a.position);
    return (diff <= positionMultiplier * k);
}
