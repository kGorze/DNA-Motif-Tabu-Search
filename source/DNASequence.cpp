//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/DNASequence.h"

DNASequence::DNASequence(const std::string &nm, const std::string &seq, const std::vector<int> &quals)
    : name(nm), bases(seq), qualityScores(quals) {}