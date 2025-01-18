#include "../include/GraphGenerator.h"
#include <random>
#include <chrono>

Graph GraphGenerator::generate(int n, double a, double b) {
    Graph G(n);
    std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(a,b);

    std::vector<double> w(n);
    for (int i = 0; i < n; i++) {
        w[i] = dist(rng);
    }

    std::uniform_real_distribution<double> edgeDist(0.0, 1.0);
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double p = (w[i] + w[j]) / 2.0;
            if (edgeDist(rng) < p) {
                G.add_edge(i, j);
            }
        }
    }

    return G;
}