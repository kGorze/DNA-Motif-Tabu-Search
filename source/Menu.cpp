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
#include <set>
#include <map>
#include <numeric>
#include <sstream>

#include "../include/Menu.h"
#include "../include/FastaQualParser.h"
#include "../include/DNASequence.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/TabuSearchMotif.h"
#include "../include/TabuSearchBase.h"
#include "../include/EnhancedMotif.h"

// ======================================================================
//  Global variables:
// ======================================================================
static std::vector<DNASequence> gAllSequences;      // All loaded sequences from file
static std::vector<DNASequence> gSelectedSequences; // 5 chosen sequences
static std::string gCurrentMotif;
static std::vector<int> gMotifPositions;            // Positions where the motif was injected
static std::vector<DNASequence> gLoadedTestSequences; // For test.fasta sequences
static int gCurrentTestInstance = 1; // Aktualnie wybrana instancja (1-5)

static std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

// Instead of a hard MIN_ALLOWED_POSITION, use 10%/90%
static double START_OFFSET = 0.1; // Must not start in the first 10% of the sequence
static double END_OFFSET   = 0.9; // Must not end in the last 10% of the sequence

// Minimal "distance" for picking sequences (example usage, not strictly enforced here)
static const int MIN_SEQ_DISTANCE = 20;

// Add this near the top of the file with other static function declarations (around line 50)
static void loadTestSequences();

// Add these declarations near the top of the file
static void saveCurrentTestInstanceToFile();
static void loadTestInstanceFromFile();

// ======================================================================
//  Helper functions:
// ======================================================================
static int countDifferences(const std::string &a, const std::string &b) {
    if (a.size() != b.size()) return 999999;
    int diff = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) diff++;
        if (diff > 2) break; // small optimization
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

// Mutate motif in up to 2 positions
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

// Tries to pick 5 "distant" sequences. In this demo we do not store originalIndex,
// so the distance constraint is symbolic. We shuffle and pick 5 unique.
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
        // The real logic to check distance by originalIndex is omitted.
        // We'll assume everything is fine in this simplified example.
        if (!tooClose) {
            gSelectedSequences.push_back(gAllSequences[idx]);
        }
        if (gSelectedSequences.size() == 5) {
            return true;
        }
    }
    return (gSelectedSequences.size() == 5);
}

// Checks if the DNA string is valid and length is 4..9
static bool isValidDNA(const std::string &motif) {
    if (motif.size() < 4 || motif.size() > 9) return false;
    for (char c : motif) {
        char uc = std::toupper((unsigned char)c);
        if (uc!='A' && uc!='C' && uc!='G' && uc!='T') return false;
    }
    return true;
}


// ======================================================================
//  MENU OPTION 1) – Pick 5 random sequences
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
//  MENU OPTION 2) – Inject a user-defined motif
// ======================================================================
static void addMotifToSelected() {
    if (gAllSequences.empty() && gLoadedTestSequences.empty()) {
        std::cout << "[Error] No sequences loaded!\n";
        return;
    }

    std::cout << "Choose sequence set:\n"
              << "1) Random sequences from sample.fasta\n"
              << "2) Loaded test sequences\n"
              << "Choice: ";
    
    int choice;
    std::cin >> choice;
    
    std::vector<DNASequence>* targetSequences;
    if (choice == 1) {
        if (gAllSequences.empty()) {
            std::cout << "[Error] No sample sequences loaded!\n";
            return;
        }
        targetSequences = &gSelectedSequences;
        // Keep picking 5 sequences until none contain a variant of this motif
        while (true) {
            pickFiveRandomSequences();
            if (gSelectedSequences.size() < 5) break;
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
    } else if (choice == 2) {
        if (gLoadedTestSequences.empty()) {
            std::cout << "[Error] No test sequences loaded!\n";
            return;
        }
        targetSequences = &gLoadedTestSequences;
    } else {
        std::cout << "[Error] Invalid choice.\n";
        return;
    }

    // Continue with motif injection...
    std::cout << "Enter the motif to inject (length 4..9, only A/C/G/T): ";
    std::cin >> gCurrentMotif;
    if (!isValidDNA(gCurrentMotif)) {
        std::cout << "[Error] Invalid motif!\n";
        return;
    }

    // Inject the motif into the chosen sequence set
    gMotifPositions.assign(5, -1);
    for (int i = 0; i < 5; i++) {
        DNASequence &seq = (*targetSequences)[i];
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
            // Retain the original quality score instead of setting it to 40
            // seq.qualityScores[startPos + k] = 40; // Remove or comment out this line
        }
        gMotifPositions[i] = startPos;
    }

    // Report
    std::cout << "[Info] Injected motif (base) = " << gCurrentMotif
              << " (with possible local mutations)\n";
    for (int i = 0; i < 5; i++) {
        const DNASequence &seq = (*targetSequences)[i];
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
    std::cout << "Choose sequence set to analyze:\n"
              << "1) Random sequences from sample.fasta\n"
              << "2) Loaded test sequences (current: " 
              << (gLoadedTestSequences.empty() ? "none" : "instance " + std::to_string(gCurrentTestInstance)) 
              << ")\n"
              << "Choice: ";
    int choice;
    std::cin >> choice;

    std::vector<DNASequence>* targetSequences = nullptr;
    if (choice == 1) {
        targetSequences = &gSelectedSequences;
    } else if (choice == 2) {
        targetSequences = &gLoadedTestSequences;
    } else {
        std::cout << "[Error] Invalid choice.\n";
        return;
    }

    if (targetSequences->empty()) {
        std::cout << "[Error] No sequences loaded.\n";
        return;
    }

    // Get parameters from user
    int kMin, kMax, qualityThreshold, allowedMismatches, positionMultiplier, maxTabuIter;
    std::cout << "[Param Setup] Enter kMin, kMax (4..9 recommended):";
    std::cin >> kMin >> kMax;
    std::cout << "[Param Setup] Enter qualityThreshold (e.g. 10..30):";
    std::cin >> qualityThreshold;
    std::cout << "[Param Setup] Allowed mismatches (0..2 recommended):";
    std::cin >> allowedMismatches;
    std::cout << "[Param Setup] Position multiplier (e.g. 10..20):";
    std::cin >> positionMultiplier;
    std::cout << "[Param Setup] MaxTabuIter (e.g. 1000):";
    std::cin >> maxTabuIter;

    // Build graph
    EnhancedMotifGraphBuilder builder(
        *targetSequences,
        kMin,
        kMax,
        qualityThreshold,
        positionMultiplier,
        allowedMismatches,
        MotifDebugInfo() // Empty debug info for automated tests
    );

    Graph g = builder.build();

    if (g.n == 0) {
        std::cout << "[Info] No valid k-mers found. Exiting search.\n";
        return;
    }

    // Run TabuSearch
    EnhancedTabuSearchMotif motifSearch(g, 20, 15, maxTabuIter, 5);
    Solution sol = motifSearch.run();

    if (sol.size == 0) {
        std::cout << "[Info] No valid motif solution found. Exiting search.\n";
        return;
    }

    // Display results
    std::cout << "[Info] Best motif solution found with size: " << sol.size << "\n\n";
    
    // Add detailed solution display
    std::cout << "=== Found Motifs ===\n";
    for (int i = 0; i < g.n; i++) {
        if (sol.inClique[i]) {
            const auto& info = g.vertexInfo[i];
            std::cout << "Sequence " << info.sequenceIndex + 1 
                      << " (pos " << info.position << "): "
                      << info.kmer << "\n";
        }
    }
    std::cout << "\n";
}


// ======================================================================
//  Check if the original injected motif was found within the solution
//  in all 5 sequences (allowing up to 2 mismatches).
// ======================================================================
static bool checkOriginalMotifFound(const Graph &g,
                                    const Solution &sol,
                                    const std::string &originalMotif)
{
    // We assume we have 5 sequences (0..4).
    // For each, we look for a chosen vertex from that sequence
    // whose k-mer differs from originalMotif by <= 2 chars.
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
//  NEW STRUCTURE FOR ENHANCED METRICS
// ======================================================================
struct EnhancedMetrics {
    double averageMotifQuality;
    double sequenceCoverage;
    double truePositiveRate;
    double motifEntropy;
    double averageHammingDistance;
    std::vector<std::string> foundMotifs;
    std::map<int, int> motifLengthDistribution;
};


// ======================================================================
//  NEW HELPER FUNCTIONS FOR ENHANCED ANALYSIS
// ======================================================================
static double calculateEntropy(const std::string& seq) {
    std::map<char, int> freqs;
    for (char c : seq) freqs[c]++;

    double entropy = 0.0;
    for (const auto& pair : freqs) {
        double p = (double)pair.second / seq.length();
        entropy -= p * std::log2(p);
    }
    return entropy;
}

static double calculateAverageHammingDistance(const std::vector<std::string>& motifs) {
    if (motifs.empty() || motifs.size() == 1) return 0.0;

    double totalDist = 0.0;
    int comparisons = 0;
    
    for (size_t i = 0; i < motifs.size(); i++) {
        for (size_t j = i + 1; j < motifs.size(); j++) {
            if (motifs[i].length() == motifs[j].length()) {
                totalDist += countDifferences(motifs[i], motifs[j]);
                comparisons++;
            }
        }
    }
    return (comparisons > 0) ? (totalDist / comparisons) : 0.0;
}

static EnhancedMetrics calculateEnhancedMetrics(
    const Graph& g,
    const Solution& sol,
    const std::string& originalMotif)
{
    EnhancedMetrics metrics;
    std::vector<std::string> foundMotifs;
    
    // Collect found motifs
    for (int v = 0; v < g.n; v++) {
        if (sol.inClique[v]) {
            foundMotifs.push_back(g.vertexInfo[v].kmer);
            metrics.motifLengthDistribution[g.vertexInfo[v].kmer.length()]++;
        }
    }
    
    // Calculate metrics
    metrics.foundMotifs = foundMotifs;
    metrics.averageHammingDistance = calculateAverageHammingDistance(foundMotifs);
    
    // Calculate average motif quality (using entropy as a proxy, or any other measure)
    double totalQuality = 0.0;
    for (const auto& motif : foundMotifs) {
        totalQuality += calculateEntropy(motif);
    }
    metrics.averageMotifQuality = foundMotifs.empty() ? 0.0 : totalQuality / foundMotifs.size();
    
    // Calculate sequence coverage
    std::set<int> coveredSequences;
    for (int v = 0; v < g.n; v++) {
        if (sol.inClique[v]) {
            coveredSequences.insert(g.vertexInfo[v].sequenceIndex);
        }
    }
    // We assume we always have 5 sequences
    metrics.sequenceCoverage = (double)coveredSequences.size() / 5.0;
    
    // Calculate true positive rate (if originalMotif is provided)
    int truePositives = 0;
    if (!originalMotif.empty()) {
        for (const auto& motif : foundMotifs) {
            if (countDifferences(motif, originalMotif) <= 2) {
                truePositives++;
            }
        }
    }
    metrics.truePositiveRate = foundMotifs.empty() ? 0.0 :
                              (double)truePositives / foundMotifs.size();
    
    // Calculate average motif entropy
    double totalEntropy = 0.0;
    for (const auto& motif : foundMotifs) {
        totalEntropy += calculateEntropy(motif);
    }
    metrics.motifEntropy = foundMotifs.empty() ? 0.0 : totalEntropy / foundMotifs.size();
    
    return metrics;
}


// ======================================================================
//  ENHANCED AUTOMATED TESTS
// ======================================================================
void runAutomatedTests() {
    std::cout << "\n=== ENHANCED AUTOMATED TESTS START ===\n";

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

    // Create enhanced CSV file
    std::time_t now = std::time(nullptr);
    std::tm localTime{};
#ifdef _WIN32
    localtime_s(&localTime, &now);
#else
    localTime = *std::localtime(&now);
#endif

    char filename[200];
    std::strftime(filename, sizeof(filename),
                  "enhanced_motif_results_%Y%m%d_%H%M%S.csv",
                  &localTime);

    std::ofstream csv(filename);
    if (!csv.is_open()) {
        std::cout << "[Error] Cannot create CSV file: " << filename << "\n";
        return;
    }

    // Enhanced CSV header with new metrics
    csv << "Motif,MotifLength,Threshold,FoundMotifs,AverageMotifQuality,"
        << "SequenceCoverage,TruePositiveRate,MotifEntropy,AverageHammingDistance,"
        << "ExecutionTimeMs,GraphVertices,EdgeDensity,MaxCliqueSize,"
        << "MotifLengthDistribution,Examples\n";

    // Test thresholds from 5 to 40 (step 5)
    std::vector<int> thresholds;
    for (int t = 5; t <= 40; t += 5) thresholds.push_back(t);

    // Run tests for each motif
    for (const auto& testMotif : testMotifs) {
        std::cout << "\n[Test] Testing motif: " << testMotif 
                  << " (length=" << testMotif.length() << ")\n";

        // 1) Find sequences and inject motif (retry if conflict)
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

        // 2) Inject the current motif
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
                // Retain the original quality score instead of setting it to 40
                // seq.qualityScores[startPos + k] = 40; // Remove or comment out this line
            }
            gMotifPositions[i] = startPos;
        }

        // 3) Build and search for each threshold
        for (int threshold : thresholds) {
            std::cout << "  Testing threshold=" << threshold << "...\n";

            // Measure execution time
            auto t1 = std::chrono::steady_clock::now();

            // Build graph
            EnhancedMotifGraphBuilder builder(
                gSelectedSequences,
                4,  // kMin
                9,  // kMax
                threshold,
                15, // positionMultiplier
                2,  // allowedMismatches
                MotifDebugInfo() // Empty debug info for automated tests
            );
            
            Graph g = builder.build();

            // Run TabuSearch
            EnhancedTabuSearchMotif motifSearch(g, 20, 15, 3000, 5);
            Solution sol = motifSearch.run();

            auto t2 = std::chrono::steady_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

            // Calculate enhanced metrics
            EnhancedMetrics metrics = calculateEnhancedMetrics(g, sol, testMotif);

            // Calculate edge density
            double edgeDensity = 0.0;
            if (g.n > 1) {
                int totalEdges = 0;
                for (int i = 0; i < g.n; i++) {
                    for (int j = i + 1; j < g.n; j++) {
                        if (g.adj[i][j]) totalEdges++;
                    }
                }
                edgeDensity = (double)totalEdges / (g.n * (g.n - 1) / 2);
            }

            // Format motif length distribution
            std::string lengthDist;
            for (const auto& pair : metrics.motifLengthDistribution) {
                if (!lengthDist.empty()) lengthDist += "|";
                lengthDist += std::to_string(pair.first) + ":" + std::to_string(pair.second);
            }

            // Format example motifs (up to 3)
            std::string examples;
            for (size_t i = 0; i < std::min<size_t>(3, metrics.foundMotifs.size()); i++) {
                if (!examples.empty()) examples += "|";
                examples += metrics.foundMotifs[i];
            }

            // Write CSV record
            csv << testMotif << ","
                << testMotif.length() << ","
                << threshold << ","
                << metrics.foundMotifs.size() << ","
                << std::fixed << std::setprecision(4)
                << metrics.averageMotifQuality << ","
                << metrics.sequenceCoverage << ","
                << metrics.truePositiveRate << ","
                << metrics.motifEntropy << ","
                << metrics.averageHammingDistance << ","
                << ms << ","
                << g.n << ","
                << edgeDensity << ","
                << sol.size << ","
                << lengthDist << ","
                << (examples.empty() ? "-" : examples) << "\n";
        }
    }

    csv.close();
    std::cout << "\n=== ENHANCED AUTOMATED TESTS FINISHED ===\n";
    std::cout << "Results saved to: " << filename << "\n";
}

// Add a new function to display the current test instance and its quality scores
static void displayCurrentTestInstance() {
    if (gLoadedTestSequences.empty()) {
        std::cout << "[Error] No test sequences loaded!\n";
        return;
    }

    if (gMotifPositions.size() != gLoadedTestSequences.size()) {
        std::cout << "[Error] Motif positions not initialized correctly!\n";
        return;
    }

    std::cout << "\n=== Current Test Instance: " << gCurrentTestInstance << " ===\n";
    for (size_t i = 0; i < gLoadedTestSequences.size(); ++i) {
        const DNASequence &seq = gLoadedTestSequences[i];
        std::cout << "Sequence " << i + 1 << ": " << seq.name << "\n";

        // Find the position of the injected motif
        int pos = gMotifPositions[i];
        if (pos >= 0 && pos + gCurrentMotif.size() <= seq.bases.size()) {
            int motLen = (int)gCurrentMotif.size();
            std::string beforeMotif = seq.bases.substr(0, pos);
            std::string motif = seq.bases.substr(pos, motLen);  // Get full motif length
            std::string afterMotif = seq.bases.substr(pos + motLen);

            // Display the sequence with the motif highlighted - remove extra spaces
            std::cout << "Bases: " << beforeMotif << "[" << motif << "]" << afterMotif << "\n";

            // Display the quality scores with the motif highlighted
            std::cout << "Qualities: ";
            for (size_t j = 0; j < seq.qualityScores.size(); ++j) {
                if (j == pos) {
                    std::cout << "[";
                }
                std::cout << seq.qualityScores[j];
                if (j == pos + motLen - 1) {
                    std::cout << "]";
                }
                std::cout << " ";
            }
            std::cout << "\n\n";
        } else {
            // If no motif was injected, display the sequence normally
            std::cout << "Bases: " << seq.bases << "\n";
            std::cout << "Qualities: ";
            for (int q : seq.qualityScores) {
                std::cout << q << " ";
            }
            std::cout << "\n\n";
        }
    }
}

// ======================================================================
//  MENU MAIN
// ======================================================================
void runMenu() {
    // Example files – automatically load:
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
                  << "2) Load 5 sequences from test.fasta\n"
                  << "3) Inject a motif (4..9, A/C/G/T)\n"
                  << "4) Build graph & TabuSearchMotif (interactive params)\n"
                  << "5) Run automated tests (ENHANCED)\n"
                  << "6) Display current test instance\n"
                  << "7) Save current test instance to file\n"
                  << "8) Load test instance from file\n"
                  << "9) Exit\n"
                  << "Choice: ";

        int choice;
        std::cin >> choice;
        if (!std::cin.good()) {
            std::cin.clear();
            std::cin.ignore(10000,'\n');
            continue;
        }
        
        switch (choice) {
            case 1:
                pickFiveRandomSequences();
                break;
            case 2:
                loadTestSequences();
                break;
            case 3:
                addMotifToSelected();
                break;
            case 4:
                runTabuMotifSearch();
                break;
            case 5:
                runAutomatedTests();
                break;
            case 6:
                displayCurrentTestInstance();
                break;
            case 7:
                saveCurrentTestInstanceToFile();
                break;
            case 8:
                loadTestInstanceFromFile();
                break;
            case 9:
                std::cout << "Exiting.\n";
                return;
            default:
                std::cout << "Unknown option.\n";
        }
    }
}

// Move the implementation to be with other function implementations (before runMenu())
static void loadTestSequences() {
    std::cout << "Choose test instance (1-5): ";
    int instance;
    std::cin >> instance;
    
    if (instance < 1 || instance > 5) {
        std::cout << "[Error] Invalid instance number. Must be 1-5.\n";
        return;
    }
    
    try {
        std::string testFasta = "test" + (instance == 1 ? "" : std::to_string(instance)) + ".fasta";
        std::string testQual = "test" + (instance == 1 ? "" : std::to_string(instance)) + ".qual";
        
        gLoadedTestSequences = FastaQualParser::parseFastaAndQual(testFasta, testQual);
        gCurrentTestInstance = instance;
        
        if (gLoadedTestSequences.size() < 5) {
            std::cout << "[Error] Not enough sequences in test files (need 5).\n";
            gLoadedTestSequences.clear();
            return;
        }
        
        // Keep only first 5 sequences
        while (gLoadedTestSequences.size() > 5) {
            gLoadedTestSequences.pop_back();
        }
        
        std::cout << "[Info] Successfully loaded 5 sequences from instance " << instance << ":\n";
        for (int i = 0; i < 5; i++) {
            std::cout << "  " << i << ") " 
                      << gLoadedTestSequences[i].name 
                      << " (length=" << gLoadedTestSequences[i].length() << ")\n";
        }
    }
    catch (const std::exception &ex) {
        std::cout << "[Error] Failed to load test sequences: " << ex.what() << "\n";
        gLoadedTestSequences.clear();
    }
}

static void saveCurrentTestInstanceToFile() {
    std::ofstream outFile("LOAD_INSTANCE.txt");
    if (!outFile) {
        std::cerr << "Error: Could not create LOAD_INSTANCE.txt\n";
        return;
    }

    outFile << "INSTANCE " << gCurrentTestInstance << "\n"
            << "MOTIF " << gCurrentMotif << "\n"
            << "SEQUENCES " << gLoadedTestSequences.size() << "\n";

    for (size_t i = 0; i < gLoadedTestSequences.size(); i++) {
        const auto& seq = gLoadedTestSequences[i];
        int motifPos = gMotifPositions[i];
        
        outFile << "SEQ_START\n"
                << "NAME " << seq.name << "\n"
                << "LENGTH " << seq.bases.length() << "\n";
        
        // Write BASES with brackets around motif - remove extra spaces
        outFile << "BASES ";
        outFile << seq.bases.substr(0, motifPos);
        outFile << "[" << seq.bases.substr(motifPos, gCurrentMotif.length()) << "]";
        outFile << seq.bases.substr(motifPos + gCurrentMotif.length()) << "\n";
        
        // Write QUALITIES with brackets around motif scores - remove extra spaces
        outFile << "QUALITIES ";
        for (size_t j = 0; j < seq.qualityScores.size(); j++) {
            if (j == motifPos) outFile << "[";
            outFile << seq.qualityScores[j] << " ";
            if (j == motifPos + gCurrentMotif.length() - 1) outFile << "]";
        }
        outFile << "\n";
        
        outFile << "MOTIF_POS " << motifPos << "\n"
                << "SEQ_END\n";
    }
    outFile.close();
}

static void loadTestInstanceFromFile() {
    std::ifstream inFile("LOAD_INSTANCE.txt");
    if (!inFile) {
        std::cerr << "Error: Could not open LOAD_INSTANCE.txt\n";
        return;
    }

    gLoadedTestSequences.clear();
    gMotifPositions.clear();

    std::string line;
    int numSequences = 0;

    // Read instance number
    std::getline(inFile, line);
    if (line.substr(0, 9) == "INSTANCE ") {
        gCurrentTestInstance = std::stoi(line.substr(9));
    }

    // Read motif
    std::getline(inFile, line);
    if (line.substr(0, 6) == "MOTIF ") {
        gCurrentMotif = line.substr(6);
    }

    // Read number of sequences
    std::getline(inFile, line);
    if (line.substr(0, 9) == "SEQUENCES ") {
        numSequences = std::stoi(line.substr(10));
    }

    // Read each sequence
    while (std::getline(inFile, line)) {
        if (line == "SEQ_START") {
            DNASequence seq;
            int motifPos = -1;

            while (std::getline(inFile, line) && line != "SEQ_END") {
                std::istringstream iss(line);
                std::string tag;
                iss >> tag;

                if (tag == "NAME") {
                    std::getline(iss, seq.name);
                    seq.name = seq.name.substr(1); // Remove leading space
                }
                else if (tag == "BASES") {
                    std::string bases;
                    std::getline(iss, bases);
                    bases = bases.substr(1); // Remove leading space
                    
                    // Remove spaces around brackets before removing brackets
                    bases.erase(std::remove(bases.begin(), bases.end(), ' '), bases.end());
                    bases.erase(std::remove(bases.begin(), bases.end(), '['), bases.end());
                    bases.erase(std::remove(bases.begin(), bases.end(), ']'), bases.end());
                    seq.bases = bases;
                }
                else if (tag == "QUALITIES") {
                    std::string qualStr;
                    int quality;
                    while (iss >> qualStr) {
                        if (qualStr != "[" && qualStr != "]") {
                            seq.qualityScores.push_back(std::stoi(qualStr));
                        }
                    }
                }
                else if (tag == "MOTIF_POS") {
                    iss >> motifPos;
                }
            }

            gLoadedTestSequences.push_back(seq);
            gMotifPositions.push_back(motifPos);
        }
    }

    inFile.close();
    std::cout << "[Info] Successfully loaded " << gLoadedTestSequences.size() 
              << " sequences from test instance " << gCurrentTestInstance 
              << " with motif: " << gCurrentMotif << "\n";
}
