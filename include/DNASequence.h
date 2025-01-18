#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <string>
#include <vector>

class DNASequence {
public:
    std::string name;                // Nazwa sekwencji (FASTA header)
    std::string bases;               // Ciąg nukleotydów
    std::vector<int> qualityScores;  // Ciąg jakości (równoległy do 'bases')

    DNASequence() = default;
    DNASequence(const std::string &nm, const std::string &seq, const std::vector<int> &quals);

    size_t length() const { return bases.size(); }
};

#endif // DNASEQUENCE_H
