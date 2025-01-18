#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

// Metadane wierzchołka dla problemu motywu
struct VertexData {
    int sequenceIndex;   // Indeks sekwencji, do której należy k-mer
    int position;        // Pozycja w sekwencji
    std::string kmer;    // Sam k-mer
};

class Graph {
public:
    int n; // liczba wierzchołków

    // Macierz sąsiedztwa: adj[u][v] = true jeśli jest krawędź (u,v)
    std::vector<std::vector<bool>> adj;

    // Dodatkowe informacje o wierzchołkach, przydatne w motif-finding
    std::vector<VertexData> vertexInfo;

    // Konstruktor - tworzy pusty graf na n wierzchołkach
    Graph(int n_);

    // Dodaje krawędź nieskierowaną (u,v)
    void add_edge(int u, int v);

    // Sprawdza, czy istnieje krawędź (u,v)
    bool is_edge(int u, int v) const;

    // Funkcja pomocnicza do budowania grafu motywu
    // Sprawdza, czy dwa wierzchołki spełniają kryteria: ten sam k-mer,
    // różne sekwencje, różnica pozycji ≤ positionMultiplier * k
    static bool canConnect(const VertexData &a, const VertexData &b, int positionMultiplier);
};

#endif // GRAPH_H
