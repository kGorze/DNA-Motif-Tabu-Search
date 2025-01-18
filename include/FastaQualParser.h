#ifndef FASTAQUALPARSER_H
#define FASTAQUALPARSER_H

#include <string>
#include <vector>
#include "DNASequence.h"

class FastaQualParser {
public:
    // Odczyt plików .fasta i .qual i zwrócenie wektora DNASequence.
    // Rzuca wyjątek std::runtime_error w razie błędów.
    static std::vector<DNASequence> parseFastaAndQual(
        const std::string &fastaFilename,
        const std::string &qualFilename
    );

private:
    static std::vector<std::string> readFasta(const std::string &fastaFilename, std::vector<std::string> &names);
    static std::vector<std::vector<int>> readQual(const std::string &qualFilename);
};

#endif // FASTAQUALPARSER_H
