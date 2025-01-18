#include "../include/DNASequence.h"

DNASequence::DNASequence(const std::string &nm, const std::string &seq, const std::vector<int> &quals)
    : name(nm), bases(seq), qualityScores(quals) {}
