#ifndef FASTAQUALPARSER_H
#define FASTAQUALPARSER_H

#include <vector>
#include <string>
#include "DNASequence.h"

class FastaQualParser {
public:
    // Wczytuje pliki FASTA i QUAL, zwraca wektor DNASequence
    // Rzuca std::runtime_error w razie problemów.
    static std::vector<DNASequence> parseFastaAndQual(
        const std::string &fastaFilename,
        const std::string &qualFilename
    );

private:
    static std::vector<std::string> readFasta(const std::string &fastaFilename,
                                              std::vector<std::string> &names);
    static std::vector<std::vector<int>> readQual(const std::string &qualFilename);
};

#endif // FASTAQUALPARSER_H
