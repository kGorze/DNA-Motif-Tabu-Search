//
// Created by konrad_guest on 18/01/2025.
//

#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <string>
#include <vector>

class DNASequence {
public:
    std::string name;                // Sequence name (from FASTA header)
    std::string bases;               // Nucleotide sequence
    std::vector<int> qualityScores;  // Quality scores (parallel to 'bases' string)

    DNASequence() = default;
    DNASequence(const std::string &nm, const std::string &seq, const std::vector<int> &quals);

    // Utility functions
    size_t length() const { return bases.size(); }
};

#endif //DNASEQUENCE_H
