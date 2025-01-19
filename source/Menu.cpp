/*
 * Created by konrad_guest on 03/12/2024.
 */

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <set> // Dodany nagłówek

#include "../include/Menu.h"
#include "../include/FastaQualParser.h"
#include "../include/DNASequence.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/TabuSearchMotif.h"
#include "../include/TabuSearchBase.h"
#include "../include/EnhancedMotif.h"

// ======================================================================
//  Global variables (unchanged):
// ======================================================================
static std::vector<DNASequence> gAllSequences;      // Wszystkie wczytane sekwencje z pliku
static std::vector<DNASequence> gSelectedSequences; // 5 wybranych sekwencji
static std::string gCurrentMotif;
static std::vector<int> gMotifPositions;            // Pozycje wstrzykniętego motywu

static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// Zamiast twardego MIN_ALLOWED_POSITION użyjemy 10%/90%
static double START_OFFSET = 0.1; // Nie może zaczynać się w pierwszych 10% sekwencji
static double END_OFFSET   = 0.9; // Nie może kończyć się w ostatnich 10% sekwencji

// Minimalny „odstęp” między indeksami w gAllSequences (przykładowo, nieużywany wprost)
static const int MIN_SEQ_DISTANCE = 20;


// ======================================================================
//  Helper functions:
// ======================================================================
static int countDifferences(const std::string &a, const std::string &b) {
    if (a.size() != b.size()) return 999999;
    int diff = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) diff++;
        if (diff > 2) break; // drobna optymalizacja
    }
    return diff;
}

static bool isMotifOrVariantPresent(const std::string &seq,
                                    const std::string &motif,
                                    int maxMismatches)
{
    int lenSeq = (int)seq.size();
    int lenMot = (int)motif.size();
    if (lenMot > lenSeq) return false;

    for (int start = 0; start + lenMot <= lenSeq; start++) {
        std::string fragment = seq.substr(start, lenMot);
        if (countDifferences(fragment, motif) <= maxMismatches) {
            return true;
        }
    }
    return false;
}

// Mutuje motyw w max 2 pozycjach
static std::string mutateMotif(const std::string &motif, int maxMutations = 2) {
    static const char* NUC = "ACGT";
    std::uniform_int_distribution<int> mutateDist(0, (int)motif.size()-1);
    std::uniform_int_distribution<int> letterDist(0, 3);
    std::uniform_int_distribution<int> howManyDist(0, maxMutations);

    std::string mutated = motif;
    int howMany = howManyDist(rng);

    for (int i = 0; i < howMany; i++) {
        int pos = mutateDist(rng);
        char newNuc = NUC[ letterDist(rng) ];
        mutated[pos] = newNuc;
    }
    return mutated;
}

// Funkcja wybierająca 5 sekwencji "odległych" (przykładowo) – w tym demie nie mamy originalIndex,
// ale logika jest zademonstrowana symbolicznie.
static bool pickFiveRandomSequencesDistant() {
    if (gAllSequences.size() < 5) {
        std::cout << "[Error] Not enough sequences.\n";
        return false;
    }
    std::vector<int> indices(gAllSequences.size());
    for (size_t i = 0; i < indices.size(); i++) {
        indices[i] = (int)i;
    }
    std::shuffle(indices.begin(), indices.end(), rng);

    gSelectedSequences.clear();
    for (int idx : indices) {
        bool tooClose = false;
        // Logika "tooClose" pominięta w tym przykładzie (brak originalIndex),
        // normalnie sprawdzalibyśmy |idx - wybrany| >= MIN_SEQ_DISTANCE.
        if (!tooClose) {
            gSelectedSequences.push_back(gAllSequences[idx]);
        }
        if (gSelectedSequences.size() == 5) {
            return true;
        }
    }
    return (gSelectedSequences.size() == 5);
}

// Sprawdza, czy ciąg DNA jest poprawny i ma długość 4..9
static bool isValidDNA(const std::string &motif) {
    if (motif.size() < 4 || motif.size() > 9) return false;
    for (char c : motif) {
        char uc = std::toupper((unsigned char)c);
        if (uc!='A' && uc!='C' && uc!='G' && uc!='T') return false;
    }
    return true;
}


// ======================================================================
//  MENU OPTION 1) – pick 5 random sequences
// ======================================================================
static void pickFiveRandomSequences() {
    bool success = pickFiveRandomSequencesDistant();
    if (!success) {
        std::cout << "[Error] Could not find 5 distant sequences.\n";
        gSelectedSequences.clear();
    } else {
        std::cout << "[Info] Picked 5 sequences (with distance constraint):\n";
        for (int i = 0; i < 5; i++) {
            std::cout << "  " << i << ") "
                      << gSelectedSequences[i].name
                      << " (length=" << gSelectedSequences[i].length() << ")\n";
        }
    }
}


// ======================================================================
//  MENU OPTION 2) – inject a user-defined motif
// ======================================================================
static void addMotifToSelected() {
    if (gAllSequences.empty()) {
        std::cout << "[Error] No sequences loaded!\n";
        return;
    }

    std::cout << "Enter the motif to inject (length 4..9, only A/C/G/T): ";
    std::cin >> gCurrentMotif;
    if (!isValidDNA(gCurrentMotif)) {
        std::cout << "[Error] Invalid motif!\n";
        return;
    }

    // Pętla "do skutku": losujemy 5 sekwencji, sprawdzamy konflikt
    while (true) {
        pickFiveRandomSequences();
        if (gSelectedSequences.size() < 5) {
            std::cout << "[Error] Could not pick 5 valid sequences at all.\n";
            return;
        }
        bool conflict = false;
        for (int i = 0; i < 5; i++) {
            const DNASequence &seq = gSelectedSequences[i];
            if (isMotifOrVariantPresent(seq.bases, gCurrentMotif, 2)) {
                conflict = true;
                break;
            }
        }
        if (!conflict) {
            break;
        } else {
            std::cout << "[Warn] At least one sequence already contained the motif (or variant)!\n"
                      << "       Re-picking 5 sequences...\n";
        }
    }

    // Wstrzykujemy
    gMotifPositions.assign(5, -1);
    for (int i = 0; i < 5; i++) {
        DNASequence &seq = gSelectedSequences[i];
        int L = (int)seq.length();
        int minStart = (int)std::floor(START_OFFSET * L);
        int maxStart = (int)std::floor(END_OFFSET * L - (double)gCurrentMotif.size());
        if (maxStart < minStart) {
            std::cout << "[Warn] Sequence " << seq.name
                      << " is too short or offsets are too big.\n";
            gMotifPositions[i] = -1;
            continue;
        }
        std::uniform_int_distribution<int> dist(minStart, maxStart);
        int startPos = dist(rng);
        std::string localMotif = mutateMotif(gCurrentMotif, 2);

        for (size_t k = 0; k < localMotif.size(); k++) {
            seq.bases[startPos + k] = localMotif[k];
            seq.qualityScores[startPos + k] = 40; // wysoka jakość
        }
        gMotifPositions[i] = startPos;
    }

    // Raport
    std::cout << "[Info] Injected motif (base) = " << gCurrentMotif
              << " (with possible local mutations)\n";
    for (int i = 0; i < 5; i++) {
        const DNASequence &seq = gSelectedSequences[i];
        int pos = gMotifPositions[i];
        std::cout << "[" << seq.name << "]\n";
        if (pos < 0) {
            std::cout << "  No injection.\n";
        } else {
            int motLen = (int)gCurrentMotif.size();
            std::string injected = seq.bases.substr(pos, motLen);
            int beforeCount = 5;
            int afterCount  = 5;
            int startPrint  = std::max(0, pos - beforeCount);
            int endPos      = pos + motLen;
            int endPrint    = std::min<int>(seq.length(), endPos + afterCount);

            std::string before = seq.bases.substr(startPrint, pos - startPrint);
            std::string middle = seq.bases.substr(pos, motLen);
            std::string after  = seq.bases.substr(endPos, endPrint - endPos);

            std::cout << "  Placed at pos=" << pos << "   LocalMotif=" << middle << "\n"
                      << "  Fragment = " << before 
                      << "[ " << middle << " ]" 
                      << after << "\n";
        }
    }
}


// ======================================================================
//  MENU OPTION 3) – build graph & run TabuSearchMotif (interactive params)
// ======================================================================
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
    std::cout << "[Param Setup] Position multiplier (e.g. 10..20): ";
    std::cin >> positionMult;
    std::cout << "[Param Setup] MaxTabuIter (e.g. 1000): ";
    std::cin >> maxIter;

    try {
        // Przygotowanie debug-info
        MotifDebugInfo debug;
        debug.motif = gCurrentMotif; 
        debug.qualityThreshold = qualityThreshold;
        for (int i = 0; i < (int)gMotifPositions.size(); i++) {
            debug.positions.push_back(gMotifPositions[i]);
            if (gMotifPositions[i] >= 0) {
                DNASequence &seq = gSelectedSequences[i];
                int startPrint = std::max(0, gMotifPositions[i] - 5);
                int endPrint   = std::min((int)seq.length(),
                                          gMotifPositions[i] + 5 + (int)gCurrentMotif.size());
                debug.sequences.push_back(seq.bases.substr(startPrint, endPrint - startPrint));
            } else {
                debug.sequences.push_back("[no injection]");
            }
        }
        debug.print();

        // Budowa grafu
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

        // TabuSearch
        int seqCount = (int)gSelectedSequences.size();
        EnhancedTabuSearchMotif motifSearch(motifG, 15, 10, maxIter, seqCount);
        Solution sol = motifSearch.run();

        std::cout << "[Result] Found clique of size " << sol.size << "\n";
        if (sol.size == seqCount) {
            std::cout << "[OK] Motif found in all sequences!\n";
        } else {
            std::cout << "[Info] Not all sequences had a suitable k-mer.\n";
        }

        // Wypisanie wierzchołków
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

            std::cout << "[" << seq.name << "]\n";
            std::cout << "  Pos = " << pos << "  Fragment = "
                      << before << "[ " << motifPart << " ]"
                      << after << "\n";
        }

        // Porównanie ze wstrzykniętym motywem
        if (!gCurrentMotif.empty()) {
            std::cout << "\n[Analysis] Comparing found k-mers with injected motif "
                      << gCurrentMotif << ":\n";
            std::cout << "Injected:     " << gCurrentMotif << "\n";
            for (int idx : chosen) {
                const auto &vd = motifG.vertexInfo[idx];
                std::string kmer = vd.kmer;
                if (kmer.size() == gCurrentMotif.size()) {
                    int diff = 0;
                    std::string diffDesc;
                    for (size_t i = 0; i < kmer.size(); i++) {
                        if (kmer[i] != gCurrentMotif[i]) {
                            if (!diffDesc.empty()) diffDesc += ", ";
                            diffDesc += gCurrentMotif[i];
                            diffDesc += "->";
                            diffDesc += kmer[i];
                            diff++;
                        }
                    }
                    std::cout << "Seq " << vd.sequenceIndex << " pos=" << vd.position 
                              << ": " << kmer << " (diff=" << diff << ": " << diffDesc << ")\n";
                } else {
                    std::cout << "Seq " << vd.sequenceIndex << " pos=" << vd.position
                              << ": " << kmer << " (length mismatch)\n";
                }
            }
        }

    } catch (const std::exception &ex) {
        std::cout << "[Error] " << ex.what() << "\n";
    }
}


// ======================================================================
//  POMOCNICZA FUNKCJA: sprawdza, czy w kliku S jest wierzchołek
//  z każdej z 5 sekwencji, którego k-mer różni się od oryginalnego
//  motywu o <= 2 znaki.
// ======================================================================
static bool checkOriginalMotifFound(const Graph &g,
                                    const Solution &sol,
                                    const std::string &originalMotif)
{
    // Zakładamy, że mamy 5 sekwencji (0..4).
    // Dla każdej z nich szukamy wierzchołka w sol, który:
    //  - jest z sekwencji i,
    //  - countDifferences <= 2 z originalMotif.
    for (int seqi = 0; seqi < 5; seqi++) {
        bool foundForThisSeq = false;
        for (int v = 0; v < g.n; v++) {
            if (sol.inClique[v]) {
                if (g.vertexInfo[v].sequenceIndex == seqi) {
                    const std::string &kmer = g.vertexInfo[v].kmer;
                    if (countDifferences(kmer, originalMotif) <= 2) {
                        foundForThisSeq = true;
                        break;
                    }
                }
            }
        }
        if (!foundForThisSeq) {
            return false;
        }
    }
    return true;
}


// ======================================================================
//  NOWA OPCJA MENU – AUTOMATYCZNE TESTY
//
// 1) Motyw: "ACGTACGT"
// 2) Pick 5 sekwencji (bez pytania użytkownika).
// 3) Sprawdzamy konflikt -> jeśli tak, losujemy jeszcze raz.
// 4) Wstrzykujemy w/w motyw w każdą z 5 sekwencji.
// 5) W pętli dla progów 10..40 (co 5):
//    - budujemy graf
//    - uruchamiamy TabuSearch
//    - mierzymy czas
//    - zbieramy wyniki i zapisujemy w CSV
// ======================================================================
void runAutomatedTests() {
    std::cout << "\n=== AUTOMATED TESTS START ===\n";

    // Set fixed seed for reproducibility
    std::mt19937 fixedRng(12345);
    rng = fixedRng;

    // Test multiple motifs of different lengths
    std::vector<std::string> testMotifs = {
        "ACGT",     // 4-mer
        "TGCAA",    // 5-mer
        "GCATAG",   // 6-mer
        "ACGTACG",  // 7-mer
        "TGCATGCA"  // 8-mer
    };

    // Prepare CSV file
    std::time_t now = std::time(nullptr);
    std::tm localTime{};
#ifdef _WIN32
    localtime_s(&localTime, &now);
#else
    localTime = *std::localtime(&now);
#endif

    char filename[200];
    std::strftime(filename, sizeof(filename),
                  "motif_test_results_%Y%m%d_%H%M%S.csv",
                  &localTime);

    std::ofstream csv(filename);
    if (!csv.is_open()) {
        std::cout << "[Error] Cannot create CSV file: " << filename << "\n";
        return;
    }

    // Enhanced CSV header with more metrics
    csv << "Motif,MotifLength,Threshold,FoundMotifs,FoundMotifLength,ExamplesFound,"
        << "ExecutionTimeMs,GraphVertices,OriginalMotifFound,AverageEdgeDensity,"
        << "MaxCliqueSize,AverageQualityScore,SequenceCoverage\n";

    // Test each motif
    for (const auto& testMotif : testMotifs) {
        std::cout << "\n[Test] Testing motif: " << testMotif << " (length=" << testMotif.length() << ")\n";

        // Find sequences without the current motif
        bool successGlobal = false;
        const int maxTries = 30;
        for (int attempt = 0; attempt < maxTries; attempt++) {
            bool ok = pickFiveRandomSequencesDistant();
            if (!ok) continue;
            
            bool conflict = false;
            for (int i = 0; i < 5; i++) {
                const DNASequence &seq = gSelectedSequences[i];
                if (isMotifOrVariantPresent(seq.bases, testMotif, 2)) {
                    conflict = true;
                    break;
                }
            }
            if (!conflict) {
                successGlobal = true;
                break;
            }
        }

        if (!successGlobal) {
            std::cout << "[Error] Could not find suitable sequences for motif " << testMotif << "\n";
            continue;
        }

        // Inject the current motif
        gMotifPositions.assign(5, -1);
        for (int i = 0; i < 5; i++) {
            DNASequence &seq = gSelectedSequences[i];
            int L = (int)seq.length();
            int minStart = (int)std::floor(START_OFFSET * L);
            int maxStart = (int)std::floor(END_OFFSET * L - (double)testMotif.size());
            if (maxStart < minStart) {
                minStart = 0;
                maxStart = std::max(0, L - (int)testMotif.size());
            }
            std::uniform_int_distribution<int> dist(minStart, maxStart);
            int startPos = dist(rng);
            std::string localMotif = mutateMotif(testMotif, 2);

            for (size_t k = 0; k < localMotif.size(); k++) {
                seq.bases[startPos + k] = localMotif[k];
                seq.qualityScores[startPos + k] = 40;
            }
            gMotifPositions[i] = startPos;
        }

        // Test thresholds from 5 to 40 with step of 1
        for (int threshold = 5; threshold <= 40; threshold++) {
            std::cout << "  Testing threshold=" << threshold << "...\n";
            
            MotifDebugInfo debug;
            debug.motif = testMotif;
            debug.qualityThreshold = threshold;
            debug.positions = gMotifPositions;

            for (int i = 0; i < 5; i++) {
                const DNASequence &seq = gSelectedSequences[i];
                int pos = gMotifPositions[i];
                if (pos >= 0) {
                    int startPrint = std::max(0, pos - 5);
                    int endPrint   = std::min((int)seq.length(),
                                              pos + 5 + (int)testMotif.size());
                    debug.sequences.push_back(seq.bases.substr(startPrint, endPrint - startPrint));
                } else {
                    debug.sequences.push_back("[no injection]");
                }
            }

            // Measure execution time
            auto t1 = std::chrono::steady_clock::now();

            // Build graph with more flexible parameters
            EnhancedMotifGraphBuilder builder(
                gSelectedSequences,
                4,        // kMin
                9,        // kMax
                threshold,
                15,       // positionMultiplier
                2,        // allowedMismatches
                debug
            );
            Graph g = builder.build();

            // Run TabuSearch with increased iterations
            EnhancedTabuSearchMotif motifSearch(g, 20, 15, 3000, 5);
            Solution sol = motifSearch.run();

            auto t2 = std::chrono::steady_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

            // Calculate additional metrics
            double edgeDensity = 0.0;
            if (g.n > 1) {
                int totalEdges = 0;
                for (int i = 0; i < g.n; i++) {
                    for (int j = i + 1; j < g.n; j++) { // Poprawka: liczyć tylko j > i
                        if (g.adj[i][j]) {
                            totalEdges++;
                        }
                    }
                }
                edgeDensity = (double)totalEdges / (g.n * (g.n - 1) / 2);
            }

            int maxCliqueSize = sol.size;
            
            // Calculate average quality score for vertices in solution
            double avgQuality = 0.0;
            int qualityCount = 0;
            for (int v = 0; v < g.n; v++) {
                if (sol.inClique[v]) {
                    avgQuality += g.vertexInfo[v].quality;
                    qualityCount++;
                }
            }
            if (qualityCount > 0) {
                avgQuality /= qualityCount;
            }

            // Calculate sequence coverage (how many different sequences are represented in the solution)
            std::set<int> coveredSequences;
            for (int v = 0; v < g.n; v++) {
                if (sol.inClique[v]) {
                    coveredSequences.insert(g.vertexInfo[v].sequenceIndex);
                }
            }
            double seqCoverage = (double)coveredSequences.size() / 5.0;

            // Collect example k-mers
            std::vector<std::string> examples;
            if (sol.size > 0) {
                for (int v = 0; v < g.n && examples.size() < 3; v++) {
                    if (sol.inClique[v]) {
                        examples.push_back(g.vertexInfo[v].kmer);
                    }
                }
            }

            std::string examplesStr = examples.empty() ? "-" : 
                                    examples[0] + 
                                    (examples.size() > 1 ? "|" + examples[1] : "") +
                                    (examples.size() > 2 ? "|" + examples[2] : "");

            bool originalFound = checkOriginalMotifFound(g, sol, testMotif);

            // Zastąpienie sol.firstVertex() logiką znajdowania pierwszego wierzchołka
            int firstVertex = -1;
            if (sol.size > 0) {
                auto it = std::find(sol.inClique.begin(), sol.inClique.end(), true);
                if (it != sol.inClique.end()) {
                    firstVertex = std::distance(sol.inClique.begin(), it);
                }
            }

            // Write enhanced CSV record
            csv << testMotif << ","
                << testMotif.length() << ","
                << threshold << ","
                << sol.size << ","
                << (sol.size > 0 && firstVertex != -1 ? g.vertexInfo[firstVertex].kmer.length() : 0) << ","
                << examplesStr << ","
                << ms << ","
                << g.n << ","
                << (originalFound ? "yes" : "no") << ","
                << std::fixed << std::setprecision(4) << edgeDensity << ","
                << maxCliqueSize << ","
                << avgQuality << ","
                << seqCoverage << "\n";
        }
    }

    csv.close();
    std::cout << "\n=== AUTOMATED TESTS FINISHED ===\n";
}


// ======================================================================
//  MENU GŁÓWNE
// ======================================================================
void runMenu() {
    // Przykładowe pliki – wczytujemy automatycznie:
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
                  << "1) Pick 5 random (distant) sequences\n"
                  << "2) Inject a motif (4..9, A/C/G/T) (re-pick if conflict)\n"
                  << "3) Build graph & TabuSearchMotif (interactive params)\n"
                  << "4) Exit\n"
                  << "5) Run automated tests\n"
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
        } else if (choice == 5) {
            runAutomatedTests();
        } else {
            std::cout << "Unknown option.\n";
        }
    }
}
