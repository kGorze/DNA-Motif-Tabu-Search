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
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

#include "../include/DNASequence.h"
#include "../include/FastaQualParser.h"
#include "../include/Graph.h"
#include "../include/GraphGenerator.h"
#include "../include/IteratedStabulus.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/Stabulus.h"
#include "../include/TabuSearchBase.h"
#include "../include/TabuSearchDeterministic.h"
#include "../include/TabuSearchMotif.h"
#include "../include/TabuSearchProbabilistic.h"

// GLOBALNE (dla uproszczenia): pula wszystkich sekwencji z duzego pliku
static std::vector<DNASequence> gAllSequences;

// 5 aktualnie wybranych sekwencji
static std::vector<DNASequence> gSelectedSequences;

// Przechowujemy parametry motywu (np. wprowadzonego sztucznie)
static std::string gCurrentMotif;
static std::vector<int> gMotifPositions; // index i -> pozycja w sekwencji i-tej

// Generator liczb losowych
static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// ------------------------------------------------------------------------------------------
// Funkcja menu 1: Losujemy 5 sekwencji z puli (gAllSequences) i zapisujemy w gSelectedSequences
// ------------------------------------------------------------------------------------------
void pickFiveRandomSequences() {
    if (gAllSequences.size() < 5) {
        std::cout << "[Blad] Za malo sekwencji w puli, wczytaj wiekszy plik!\n";
        return;
    }

    // Losujemy 5 unikalnych indeksow
    std::vector<int> indices(gAllSequences.size());
    for (size_t i = 0; i < indices.size(); i++) {
        indices[i] = (int)i;
    }
    std::shuffle(indices.begin(), indices.end(), rng);

    gSelectedSequences.clear();
    for (int i = 0; i < 5; i++) {
        gSelectedSequences.push_back(gAllSequences[indices[i]]);
    }

    // Zerujemy stare informacje o wstrzyknietym motywie
    gCurrentMotif.clear();
    gMotifPositions.clear();

    std::cout << "Wylosowano 5 sekwencji:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << "  " << i << ") " << gSelectedSequences[i].name 
                  << " (len=" << gSelectedSequences[i].length() << ")\n";
    }
}

// ------------------------------------------------------------------------------------------
// Funkcja menu 2: Dodajemy (wstrzykujemy) wskazany motyw do kazdej z 5 sekwencji
// w losowej (lub pytanej) pozycji. Wypisujemy w formacie:
// {nazwa}{pozycja}{...kilka poprzednich...}..{motyw}..{kilka nastepnych}
// ------------------------------------------------------------------------------------------
void addMotifToSelected() {
    if (gSelectedSequences.size() < 5) {
        std::cout << "[Blad] Nie wybrano 5 sekwencji.\n";
        return;
    }

    std::cout << "Podaj motyw, ktory chcesz wstrzyknac (np. 'ACGTACGT'):\n> ";
    std::cin >> gCurrentMotif;
    if (gCurrentMotif.empty()) {
        std::cout << "[Blad] Motyw jest pusty!\n";
        return;
    }

    gMotifPositions.clear();
    gMotifPositions.resize(5, -1);

    // Dla kazdej z 5 sekwencji wstrzykujemy motyw w losowym miejscu
    // (o ile motyw sie zmiesci)
    for (int i = 0; i < 5; i++) {
        DNASequence &seq = gSelectedSequences[i];
        int maxStart = (int)seq.length() - (int)gCurrentMotif.size();
        if (maxStart < 0) {
            // Jesli sekwencja krotsza od motywu, pomijamy
            std::cout << "[Uwaga] Sekwencja " << seq.name << " za krotka, nie wstrzyknieto.\n";
            gMotifPositions[i] = -1;
            continue;
        }
        std::uniform_int_distribution<int> dist(0, maxStart);
        int startPos = dist(rng);

        // Wstrzykniecie: zamieniamy bazy w zakresie [startPos, startPos+motif.size())
        for (size_t k = 0; k < gCurrentMotif.size(); k++) {
            seq.bases[startPos + k] = gCurrentMotif[k];
            // Na potrzeby testu: ustawmy qualityScores na wysokie
            seq.qualityScores[startPos + k] = 35; 
        }
        gMotifPositions[i] = startPos;
    }

    // Teraz wypisujemy w żądanym formacie
    std::cout << "Wstrzykniety motyw: " << gCurrentMotif << "\n\n";
    for (int i = 0; i < 5; i++) {
        const DNASequence &seq = gSelectedSequences[i];
        int pos = gMotifPositions[i];
        std::cout << "[" << seq.name << "]\n";

        if (pos < 0) {
            std::cout << "  Motyw nie wstrzykniety (za krotka sekwencja).\n";
            continue;
        }

        // Kilka nukleotydów przed
        int beforeCount = 5;
        int startPrint = std::max(0, pos - beforeCount);
        std::string before = seq.bases.substr(startPrint, pos - startPrint);

        // Sam motyw
        std::string motifPart = seq.bases.substr(pos, gCurrentMotif.size());

        // Kilka nukleotydów po
        int afterCount = 5;
        int endPos = pos + (int)gCurrentMotif.size();
        int endPrint = std::min((int)seq.length(), endPos + afterCount);
        std::string after = seq.bases.substr(endPos, endPrint - endPos);

        std::cout << "  Pozycja = " << pos << "\n";
        std::cout << "  Fragment = " << before << ".." << motifPart << ".." << after << "\n";
    }
}

// ------------------------------------------------------------------------------------------
// Funkcja menu 3: Budujemy graf motywu na podstawie gSelectedSequences
// i uruchamiamy TabuSearchMotif, a nastepnie analizujemy wyniki
// ------------------------------------------------------------------------------------------
void runTabuMotifSearch() {
    if (gSelectedSequences.size() < 5) {
        std::cout << "[Blad] Nie mamy 5 sekwencji do analizy.\n";
        return;
    }

    // Przykładowe parametry
    int kMin = 4;
    int kMax = 9;
    int qualityThreshold = 14;  // np. z raportu: testowanie roznych wartosci: 14, 19, 24 itp.
    int positionMultiplier = 10;

    // Budujemy graf
    MotifGraphBuilder builder(gSelectedSequences, kMin, kMax, qualityThreshold, positionMultiplier);
    Graph motifGraph = builder.build();
    int seqCount = (int)gSelectedSequences.size();

    std::cout << "Zbudowano graf motywu. Liczba wierzcholkow = " << motifGraph.n << "\n";

    // Odpalamy TabuSearchMotif - chcemy klike rozmiaru = 5
    TabuSearchMotif motifSearch(motifGraph, 10, 5, 1000, seqCount);
    Solution sol = motifSearch.run();

    // Wyświetlamy wyniki
    std::cout << "TabuSearchMotif: uzyskana klika rozmiaru " << sol.size << ".\n";
    if (sol.size == seqCount) {
        std::cout << "-> Znaleziono idealny motyw we wszystkich sekwencjach.\n";
    } else {
        std::cout << "-> Nie jest to kompletna reprezentacja we wszystkich 5 sekwencjach.\n";
    }

    std::vector<int> chosenVerts;
    for (int i = 0; i < motifGraph.n; i++) {
        if (sol.inClique[i]) chosenVerts.push_back(i);
    }

    // Prezentacja, które k-mery (sekwencja, pozycja, sam k-mer) zostały znalezione
    std::cout << "Znaleziona klika:\n";
    for (int idx : chosenVerts) {
        const VertexData &vd = motifGraph.vertexInfo[idx];
        std::cout << "  SeqIndex=" << vd.sequenceIndex
                  << " Position=" << vd.position
                  << " K-mer=" << vd.kmer << "\n";
    }

    // Dodatkowo sprawdzimy, czy któraś z pozycji pokrywa się z wstrzykniętym motywem
    if (!gCurrentMotif.empty()) {
        std::cout << "\nAnaliza sztucznego motywu: " << gCurrentMotif << "\n";
        for (int idx : chosenVerts) {
            const VertexData &vd = motifGraph.vertexInfo[idx];
            // Sprawdz, czy to "sztuczne" wystąpienie
            // tzn. czy wstrzyknęliśmy i czy te pozycje się pokrywają
            if (vd.kmer == gCurrentMotif) {
                // weryfikacja, czy position == gMotifPositions[seqIndex], zakładając dokładną zgodność
                if (vd.position == gMotifPositions[vd.sequenceIndex]) {
                    std::cout << "  -> K-mer w sekwencji " << vd.sequenceIndex
                              << " to SZTUCZNIE wprowadzony motyw!\n";
                } else {
                    std::cout << "  -> K-mer w sekwencji " << vd.sequenceIndex
                              << " to raczej NATURALNE wystąpienie " << vd.kmer << "\n";
                }
            }
        }
    }
}

// ------------------------------------------------------------------------------------------
// Główna funkcja z pętlą menu
// ------------------------------------------------------------------------------------------
int main(int argc, char** argv) {
    std::string fastaFile = "sample.fasta";
    std::string qualFile = "sample.qual";

    try {
        // Wczytanie plików sample z bieżącego katalogu
        gAllSequences = FastaQualParser::parseFastaAndQual(fastaFile, qualFile);

        std::cout << "Wczytano " << gAllSequences.size() << " sekwencji z plikow " 
                  << fastaFile << " i " << qualFile << ".\n";
    }
    catch(const std::exception &ex) {
        std::cerr << "[BLAD] Nie mozna wczytac plikow " << fastaFile << " lub " << qualFile 
                  << ": " << ex.what() << "\n";
        return 1;
    }

    // Proste menu
    while (true) {
        std::cout << "\n=== MENU ===\n"
                  << "1) Wylosuj 5 sekwencji z puli\n"
                  << "2) Dodaj motyw do wybranych sekwencji i wyswietl\n"
                  << "3) Wyszukaj motyw TabuSearch (1 wierzcholek/sekwencja)\n"
                  << "4) [Wyjscie]\n"
                  << "Wybierz opcje: ";

        int opt;
        std::cin >> opt;
        if (!std::cin.good()) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }

        switch(opt) {
        case 1:
            pickFiveRandomSequences();
            break;
        case 2:
            addMotifToSelected();
            break;
        case 3:
            runTabuMotifSearch();
            break;
        case 4:
            std::cout << "Koniec.\n";
            return 0;
        default:
            std::cout << "Nieznana opcja.\n";
            break;
        }
    }

    return 0;
}


