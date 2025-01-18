#include "../include/FastaQualParser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

std::vector<DNASequence> FastaQualParser::parseFastaAndQual(
    const std::string &fastaFilename,
    const std::string &qualFilename
) {
    // Wczytaj FASTA
    std::vector<std::string> names;
    std::vector<std::string> sequences = readFasta(fastaFilename, names);

    // Wczytaj QUAL
    std::vector<std::vector<int>> quals = readQual(qualFilename);

    if (sequences.size() != quals.size()) {
        throw std::runtime_error("Liczba sekwencji w FASTA != liczba w QUAL!");
    }

    // Łączymy w DNASequence
    std::vector<DNASequence> result;
    result.reserve(sequences.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].size() != quals[i].size()) {
            throw std::runtime_error("Długość sekwencji != liczba jakości w sekwencji nr " + std::to_string(i));
        }
        result.emplace_back(names[i], sequences[i], quals[i]);
    }
    return result;
}

std::vector<std::string> FastaQualParser::readFasta(const std::string &fastaFilename,
                                                    std::vector<std::string> &names) {
    std::ifstream in(fastaFilename);
    if (!in.is_open()) {
        throw std::runtime_error("Nie mogę otworzyć pliku FASTA: " + fastaFilename);
    }

    std::vector<std::string> seqs;
    std::string line, currentSeq, currentName;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentSeq.empty()) {
                seqs.push_back(currentSeq);
                names.push_back(currentName);
            }
            currentName = line.substr(1); // usuń '>'
            currentSeq.clear();
        } else {
            // fragment sekwencji
            currentSeq += line;
        }
    }
    if (!currentSeq.empty()) {
        seqs.push_back(currentSeq);
        names.push_back(currentName);
    }
    return seqs;
}

std::vector<std::vector<int>> FastaQualParser::readQual(const std::string &qualFilename) {
    std::ifstream in(qualFilename);
    if (!in.is_open()) {
        throw std::runtime_error("Nie mogę otworzyć pliku QUAL: " + qualFilename);
    }

    std::vector<std::vector<int>> allQuals;
    std::string line;
    std::vector<int> currentQuals;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentQuals.empty()) {
                allQuals.push_back(currentQuals);
            }
            currentQuals.clear();
        } else {
            // odczytujemy liczby
            std::stringstream ss(line);
            int q;
            while (ss >> q) {
                currentQuals.push_back(q);
            }
        }
    }
    if (!currentQuals.empty()) {
        allQuals.push_back(currentQuals);
    }
    return allQuals;
}
