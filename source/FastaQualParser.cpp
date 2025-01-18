//
// Created by konrad_guest on 18/01/2025.
//
#include "../include/FastaQualParser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

std::vector<DNASequence> FastaQualParser::parseFastaAndQual(
    const std::string &fastaFilename,
    const std::string &qualFilename
) {
    // Read FASTA
    std::vector<std::string> names;
    std::vector<std::string> sequences = readFasta(fastaFilename, names);

    // Read QUAL
    std::vector<std::vector<int>> qualities = readQual(qualFilename);

    if (sequences.size() != qualities.size()) {
        throw std::runtime_error("Mismatch between number of FASTA sequences and QUAL sequences.");
    }

    // Combine into DNASequence objects
    std::vector<DNASequence> result;
    result.reserve(sequences.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].size() != qualities[i].size()) {
            throw std::runtime_error("Sequence length does not match number of quality scores for sequence index " + std::to_string(i));
        }
        result.emplace_back(names[i], sequences[i], qualities[i]);
    }
    return result;
}

std::vector<std::string> FastaQualParser::readFasta(const std::string &fastaFilename, std::vector<std::string> &names) {
    std::ifstream in(fastaFilename);
    if (!in.good()) {
        throw std::runtime_error("Cannot open FASTA file: " + fastaFilename);
    }

    std::vector<std::string> sequences;
    std::string line, currentSeq;
    std::string currentName;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // new sequence header
            if (!currentSeq.empty()) {
                sequences.push_back(currentSeq);
                names.push_back(currentName);
            }
            currentName = line.substr(1); // strip '>'
            currentSeq.clear();
        } else {
            // part of sequence
            currentSeq += line;
        }
    }
    // last sequence
    if (!currentSeq.empty()) {
        sequences.push_back(currentSeq);
        names.push_back(currentName);
    }

    return sequences;
}

std::vector<std::vector<int>> FastaQualParser::readQual(const std::string &qualFilename) {
    std::ifstream in(qualFilename);
    if (!in.good()) {
        throw std::runtime_error("Cannot open QUAL file: " + qualFilename);
    }

    std::vector<std::vector<int>> allQuals;
    std::string line;
    std::vector<int> currentQuals;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // new sequence
            if (!currentQuals.empty()) {
                allQuals.push_back(currentQuals);
            }
            currentQuals.clear();
        } else {
            // parse quality scores
            std::stringstream ss(line);
            int score;
            while (ss >> score) {
                currentQuals.push_back(score);
            }
        }
    }

    // last sequence
    if (!currentQuals.empty()) {
        allQuals.push_back(currentQuals);
    }

    return allQuals;
}
