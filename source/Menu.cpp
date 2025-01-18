#include "../include/Menu.h"
#include "../include/FastaQualParser.h"
#include "../include/DNASequence.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/TabuSearchMotif.h"
#include "../include/TabuSearchBase.h"

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <string>

// ======================== GLOBALNE W TYM PLIKU ========================
static std::vector<DNASequence> gAllSequences;      // Wszystkie wczytane sekwencje
static std::vector<DNASequence> gSelectedSequences; // 5 wybranych sekwencji
static std::string gCurrentMotif;
static std::vector<int> gMotifPositions;            // Pozycje wstrzykniętego motywu

static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// Minimalny dozwolony start (nie chcemy, by motyw zaczynał się np. tuż przy 0).
static const int MIN_ALLOWED_POSITION = 5;

static bool isValidDNA(const std::string &motif) {
    if (motif.size() < 4 || motif.size() > 9) return false;
    for (char c : motif) {
        char uc = std::toupper((unsigned char)c);
        if (uc!='A' && uc!='C' && uc!='G' && uc!='T') return false;
    }
    return true;
}

// ----------------------------------------------------
// 1) Losowy wybór 5 sekwencji z puli
// ----------------------------------------------------
static void pickFiveRandomSequences() {
    if (gAllSequences.size() < 5) {
        std::cout << "[Error] Not enough sequences to pick 5.\n";
        return;
    }
    std::vector<int> indices(gAllSequences.size());
    for (size_t i = 0; i < indices.size(); i++) {
        indices[i] = (int)i;
    }
    std::shuffle(indices.begin(), indices.end(), rng);

    gSelectedSequences.clear();
    for (int i = 0; i < 5; i++) {
        gSelectedSequences.push_back(gAllSequences[ indices[i] ]);
    }
    gCurrentMotif.clear();
    gMotifPositions.assign(5, -1);

    std::cout << "[Info] Picked 5 random sequences:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << "  " << i << ") " 
                  << gSelectedSequences[i].name
                  << " (length=" << gSelectedSequences[i].length() << ")\n";
    }
}

// ----------------------------------------------------
// 2) Wstrzyknięcie motywu do każdej z 5 sekwencji
// ----------------------------------------------------
static void addMotifToSelected() {
    if (gSelectedSequences.size() < 5) {
        std::cout << "[Error] Please pick 5 sequences first (option 1)!\n";
        return;
    }

    std::cout << "Enter the motif to inject (length 4..9, only A/C/G/T): ";
    std::cin >> gCurrentMotif;
    if (!isValidDNA(gCurrentMotif)) {
        std::cout << "[Error] Invalid motif!\n";
        return;
    }

    gMotifPositions.assign(5, -1);

    for (int i = 0; i < 5; i++) {
        DNASequence &seq = gSelectedSequences[i];
        int maxStart = (int)seq.length() - (int)gCurrentMotif.size();
        if (maxStart < MIN_ALLOWED_POSITION) {
            std::cout << "[Warn] Seq " << seq.name 
                      << " is too short to place motif at >= " 
                      << MIN_ALLOWED_POSITION << "\n";
            gMotifPositions[i] = -1;
            continue;
        }

        // Wybieramy losową pozycję >= MIN_ALLOWED_POSITION
        static const int MAX_ATTEMPTS = 1000;
        std::uniform_int_distribution<int> dist(0, maxStart);

        int startPos = -1;
        for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
            int cand = dist(rng);
            if (cand >= MIN_ALLOWED_POSITION) {
                startPos = cand;
                break;
            }
        }
        if (startPos < 0) {
            std::cout << "[Warn] Could not place motif in seq " 
                      << seq.name << "\n";
            continue;
        }

        // Nadpisujemy fragment bazy
        for (size_t k = 0; k < gCurrentMotif.size(); k++) {
            seq.bases[startPos + k] = gCurrentMotif[k];
            // Ustawiamy wysoką jakość, aby nie zostało wycięte
            seq.qualityScores[startPos + k] = 40;
        }
        gMotifPositions[i] = startPos;
    }

    // Wypisujemy wyniki
    std::cout << "[Info] Injected motif: " << gCurrentMotif << "\n";
    for (int i = 0; i < 5; i++) {
        const DNASequence &seq = gSelectedSequences[i];
        int pos = gMotifPositions[i];
        std::cout << "[" << seq.name << "]\n";
        if (pos < 0) {
            std::cout << "  No injection.\n";
        } else {
            std::cout << "  Placed at pos=" << pos << "\n";
            int beforeCount = 5;
            int startPrint = std::max(0, pos - beforeCount);
            std::string before = seq.bases.substr(startPrint, pos - startPrint);

            std::string motifPart = seq.bases.substr(pos, gCurrentMotif.size());

            int afterCount = 5;
            int endPos = pos + (int)gCurrentMotif.size();
            int endPrint = std::min<int>((int)seq.length(), endPos + afterCount);
            std::string after = seq.bases.substr(endPos, endPrint - endPos);

            std::cout << "  Fragment = " << before 
                      << ".." << motifPart
                      << ".." << after << "\n";
        }
    }
}

// ----------------------------------------------------
// 3) Budowa grafu & uruchomienie TabuSearchMotif
//    Dodana możliwość wprowadzenia własnych parametrów
// ----------------------------------------------------
static void runTabuMotifSearch() {
    if (gSelectedSequences.size() < 5) {
        std::cout << "[Error] No 5 sequences selected yet.\n";
        return;
    }

    int kMin, kMax, qualityThreshold, mismatch, positionMult, maxIter;
    std::cout << "\n[Param Setup] Enter kMin, kMax (4..9 recommended): ";
    std::cin >> kMin >> kMax;
    std::cout << "[Param Setup] Enter qualityThreshold (e.g. 10..30): ";
    std::cin >> qualityThreshold;
    std::cout << "[Param Setup] Allowed mismatches (0..2 recommended): ";
    std::cin >> mismatch;
    std::cout << "[Param Setup] Position multiplier (np. 10..20): ";
    std::cin >> positionMult;
    std::cout << "[Param Setup] MaxTabuIter (np. 1000): ";
    std::cin >> maxIter;

    std::cout << "[Debug] Building graph with kMin=" << kMin 
              << ", kMax=" << kMax 
              << ", mismatch=" << mismatch
              << ", qualityThr=" << qualityThreshold
              << ", posMult=" << positionMult << "\n";

    try {
        // ---------------------------------------------------
        // PRZYGOTOWANIE STRUKTURY DEBUGUJĄCEJ
        // ---------------------------------------------------
        MotifDebugInfo debug;
        debug.motif = gCurrentMotif;  // np. z ostatniego wstrzyknięcia
        debug.qualityThreshold = qualityThreshold;

        // Pozycje i fragmenty z wstrzyknięcia (o ile przechowujemy)
        // gMotifPositions to wektor o długości 5 (dla 5 sekwencji),
        // w którym są zapisywane pozycje wstrzyknięcia.
        for (int i = 0; i < (int)gMotifPositions.size(); i++) {
            debug.positions.push_back(gMotifPositions[i]);
            // Zapiszmy też fragment sekwencji z danej pozycji
            if (gMotifPositions[i] >= 0) {
                DNASequence &seq = gSelectedSequences[i];
                // Zapisz kawałek sekwencji do debugowania
                int startPrint = std::max(0, gMotifPositions[i] - 5);
                int endPrint   = std::min((int)seq.length(), gMotifPositions[i] + 5 + (int)gCurrentMotif.size());
                debug.sequences.push_back(seq.bases.substr(startPrint, endPrint - startPrint));
            } else {
                debug.sequences.push_back("[no injection]");
            }
        }

        // Można wydrukować info debugowe w tym miejscu
        debug.print();

        // ---------------------------------------------------
        // UŻYCIE ULEPSZONEGO BUDOWNICZEGO GRAFU
        // ---------------------------------------------------
        EnhancedMotifGraphBuilder builder(
            gSelectedSequences,
            kMin,
            kMax,
            qualityThreshold,
            positionMult,
            mismatch,
            debug
        );

        Graph motifG = builder.build();
        std::cout << "[Debug] Graph has " << motifG.n << " vertices.\n";

        // ---------------------------------------------------
        // URUCHOMIENIE TABU SEARCH
        // ---------------------------------------------------
        int seqCount = (int)gSelectedSequences.size();
        // Zamiast TabuSearchMotif używamy EnhancedTabuSearchMotif
        EnhancedTabuSearchMotif motifSearch(motifG, 15, 10, maxIter, seqCount);
        Solution sol = motifSearch.run();

        std::cout << "[Result] Found clique of size " << sol.size << "\n";
        if (sol.size == seqCount) {
            std::cout << "[OK] Motif found in all sequences!\n";
        } else {
            std::cout << "[Info] Not all sequences had a suitable k-mer.\n";
        }

        // Wypisujemy wierzchołki kliki z wyróżnieniem i kontekstem
        std::vector<int> chosen;
        for (int i = 0; i < motifG.n; i++) {
            if (sol.inClique[i]) {
                chosen.push_back(i);
            }
        }
        for (int idx : chosen) {
            const auto &vd = motifG.vertexInfo[idx];
            const DNASequence &seq = gSelectedSequences[vd.sequenceIndex];
            int pos = vd.position;
            int kmerLength = (int)vd.kmer.size();

            int beforeCount = 5;
            int afterCount  = 5;

            int startPrint = std::max(0, pos - beforeCount);
            std::string before = seq.bases.substr(startPrint, pos - startPrint);

            std::string motifPart = seq.bases.substr(pos, kmerLength);

            int endPos   = pos + kmerLength;
            int endPrint = std::min((int)seq.length(), endPos + afterCount);
            std::string after = seq.bases.substr(endPos, endPrint - endPos);

            // Wyświetlenie wyników w nowym formacie z wyróżnieniem motywu i indeksu
            std::cout << "[" << seq.name << "]\n";
            std::cout << "  Pos = " << pos << "  Fragment = "
                      << before << "[ " << motifPart << " ]"
                      << after << "\n";
        }

        // Sprawdzanie, czy znaleziony k-mer jest zbliżony do wstrzykniętego motywu
        // Replace the existing analysis code at the end of runTabuMotifSearch() with:

        // Sprawdzanie, czy znaleziony k-mer jest zbliżony do wstrzykniętego motywu
        if (!gCurrentMotif.empty()) {
            std::cout << "\n[Analysis] Comparing found k-mers with injected motif " 
                      << gCurrentMotif << ":\n";
    
            // First print the injected motif as reference
            std::cout << "Injected:     " << gCurrentMotif << "\n";
    
            // Store all kmers for alignment
            std::vector<std::string> foundKmers;
    
            // Print comparison for each sequence
            for (int idx : chosen) {
                const auto &vd = motifG.vertexInfo[idx];
                std::string kmer = vd.kmer;
                foundKmers.push_back(kmer);
        
                // Only compare if lengths match
                if (kmer.size() == gCurrentMotif.size()) {
                    std::cout << "Seq " << vd.sequenceIndex 
                              << " (" << vd.position << "):   "
                              << kmer << "  (diff=";
            
                    // Count and collect differences
                    int diff = 0;
                    std::string diffDesc;
                    for (size_t i = 0; i < gCurrentMotif.size(); i++) {
                        if (kmer[i] != gCurrentMotif[i]) {
                            if (!diffDesc.empty()) diffDesc += ", ";
                            diffDesc += gCurrentMotif[i];
                            diffDesc += "->";
                            diffDesc += kmer[i];
                            diff++;
                        }
                    }
                    std::cout << diff << ": " << diffDesc << ")\n";
                }
            }
    
            // Print alignment of found k-mers
            std::cout << "\nAlignment of found k-mers:\n";
            for (const auto& kmer : foundKmers) {
                std::cout << kmer << "\n";
            }
        }
    }
    catch (const std::exception &ex) {
        std::cout << "[Error] " << ex.what() << "\n";
    }
}

// ----------------------------------------------------
// 4) Menu główne
// ----------------------------------------------------
void runMenu() {
    // Wczytujemy większy plik "sample.fasta" i "sample.qual" 
    // i wyciągamy z niego sekwencje. Następnie wybieramy z nich 5.
    std::string fastaFile = "sample.fasta";
    std::string qualFile  = "sample.qual";

    try {
        gAllSequences = FastaQualParser::parseFastaAndQual(fastaFile, qualFile);
        std::cout << "[Info] Loaded " << gAllSequences.size() 
                  << " sequences from " << fastaFile 
                  << " + " << qualFile << "\n";
    }
    catch (const std::exception &ex) {
        std::cout << "[Error] " << ex.what() << "\n";
        return;
    }

    while (true) {
        std::cout << "\n=== MENU ===\n"
                  << "1) Pick 5 random sequences\n"
                  << "2) Inject a motif (4..9, A/C/G/T)\n"
                  << "3) Build graph & TabuSearchMotif (interactive params)\n"
                  << "4) Exit\n"
                  << "Choice: ";
        int choice;
        std::cin >> choice;
        if (!std::cin.good()) {
            std::cin.clear();
            std::cin.ignore(10000,'\n');
            continue;
        }
        if (choice == 1) {
            pickFiveRandomSequences();
        } else if (choice == 2) {
            addMotifToSelected();
        } else if (choice == 3) {
            runTabuMotifSearch();
        } else if (choice == 4) {
            std::cout << "Exiting.\n";
            return;
        } else {
            std::cout << "Unknown option.\n";
        }
    }
}
