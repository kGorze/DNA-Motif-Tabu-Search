#include "../include/Graph.h"

Graph::Graph(int n_) : n(n_), adj(n_, std::vector<bool>(n_, false)), vertexInfo(n_) {}

void Graph::add_edge(int u, int v) {
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
    // 1) k-mery muszą być identyczne
    // 2) Różne sekwencje
    // 3) Różnica pozycji ≤ positionMultiplier * długość k-meru
    if (a.sequenceIndex == b.sequenceIndex) return false;
    if (a.kmer != b.kmer) return false;
    int k = (int)a.kmer.size();
    int diff = (a.position > b.position)
               ? (a.position - b.position)
               : (b.position - a.position);
    return (diff <= positionMultiplier * k);
}
