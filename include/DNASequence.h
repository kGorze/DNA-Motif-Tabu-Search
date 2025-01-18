#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <string>
#include <vector>

class DNASequence {
public:
    std::string name;                // Nazwa sekwencji (z nagłówka FASTA)
    std::string bases;               // Sekwencja nukleotydów
    std::vector<int> qualityScores;  // Jakości (równoległe do 'bases')

    DNASequence() = default;
    DNASequence(const std::string &nm, const std::string &seq, const std::vector<int> &quals);

    // Długość sekwencji
    size_t length() const { return bases.size(); }
};

#endif // DNASEQUENCE_H
