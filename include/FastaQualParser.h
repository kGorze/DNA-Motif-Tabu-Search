//
// Created by konrad_guest on 18/01/2025.
//

#ifndef FASTAQUALPARSER_H
#define FASTAQUALPARSER_H

#include <string>
#include <vector>
#include "DNASequence.h"

class FastaQualParser {
public:
    // Reads a FASTA file and corresponding QUAL file
    // Returns a vector of DNASequence objects
    // If the .qual file is missing or not matched, throws std::runtime_error
    static std::vector<DNASequence> parseFastaAndQual(
        const std::string &fastaFilename,
        const std::string &qualFilename
    );

private:
    static std::vector<std::string> readFasta(const std::string &fastaFilename, std::vector<std::string> &names);
    static std::vector<std::vector<int>> readQual(const std::string &qualFilename);
};


#endif //FASTAQUALPARSER_H
