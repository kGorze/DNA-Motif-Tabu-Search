#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

struct VertexData {
    int sequenceIndex;  // Numer sekwencji
    int position;       // Pozycja w oryginalnej sekwencji
    std::string kmer;   // Tekst k-mera
    int quality;        // Średnia jakość k-meru
};

class Graph {
public:
    int n; // liczba wierzchołków
    std::vector<std::vector<bool>> adj;
    std::vector<VertexData> vertexInfo;

    Graph(int n_);

    void add_edge(int u, int v);
    bool is_edge(int u, int v) const;

    // canConnect z uwzględnieniem limitu mismatch, etc.
    static bool canConnect(const VertexData &a,
                           const VertexData &b,
                           int positionMultiplier,
                           int allowedMismatches);

private:
    // Pomocnicza funkcja do porównywania k-merów z dozwoloną liczbą różnic
    static bool isWithinMismatch(const std::string &s1,
                                 const std::string &s2,
                                 int maxMismatch);
};

#endif // GRAPH_H