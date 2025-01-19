#include "../include/Graph.h"
#include <cmath>
#include <stdexcept>

static int countDifferences(const std::string &s1, const std::string &s2) {
    if (s1.size() != s2.size()) return 999999; // large => fail
    int diff = 0;
    for (size_t i = 0; i < s1.size(); i++) {
        if (s1[i] != s2[i]) diff++;
    }
    return diff;
}

Graph::Graph(int n_)
: n(n_),
  adj(n_, std::vector<bool>(n_, false)),
  vertexInfo(n_) 
{}

void Graph::add_edge(int u, int v) {
    if (u < 0 || v < 0 || u >= n || v >= n) return;
    if (u == v) return;

    adj[u][v] = true;
    adj[v][u] = true;
}

bool Graph::is_edge(int u, int v) const {
    if (u < 0 || v < 0 || u >= n || v >= n) return false;
    return adj[u][v];
}

bool Graph::canConnect(const VertexData &a,
                       const VertexData &b,
                       int positionMultiplier,
                       int allowedMismatches)
{
    if (a.sequenceIndex == b.sequenceIndex) {
        return false;
    }

    int diff = countDifferences(a.kmer, b.kmer);
    if (diff > allowedMismatches) {
        return false;
    }

    int k = (int)a.kmer.size();
    int dPos = std::abs(a.position - b.position);
    if (dPos > positionMultiplier * k) {
        return false;
    }

    return true;
}