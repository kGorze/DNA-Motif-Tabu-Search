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
    std::vector<std::vector<int>> qualities = readQual(qualFilename);

    if (sequences.size() != qualities.size()) {
        throw std::runtime_error("Liczba sekwencji w FASTA nie odpowiada liczbie sekwencji w QUAL!");
    }

    // Połącz w DNASequence
    std::vector<DNASequence> result;
    result.reserve(sequences.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].size() != qualities[i].size()) {
            throw std::runtime_error(
                "Długość sekwencji nie pokrywa się z liczbą wartości QUAL dla sekwencji nr " + std::to_string(i)
            );
        }
        result.emplace_back(names[i], sequences[i], qualities[i]);
    }
    return result;
}

std::vector<std::string> FastaQualParser::readFasta(const std::string &fastaFilename, std::vector<std::string> &names) {
    std::ifstream in(fastaFilename);
    if (!in.good()) {
        throw std::runtime_error("Nie można otworzyć pliku FASTA: " + fastaFilename);
    }

    std::vector<std::string> sequences;
    std::string line, currentSeq;
    std::string currentName;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // Nagłówek nowej sekwencji
            if (!currentSeq.empty()) {
                sequences.push_back(currentSeq);
                names.push_back(currentName);
            }
            currentName = line.substr(1); // usunięcie '>'
            currentSeq.clear();
        } else {
            // Kolejne fragmenty sekwencji
            currentSeq += line;
        }
    }
    // Ostatnia sekwencja
    if (!currentSeq.empty()) {
        sequences.push_back(currentSeq);
        names.push_back(currentName);
    }

    return sequences;
}

std::vector<std::vector<int>> FastaQualParser::readQual(const std::string &qualFilename) {
    std::ifstream in(qualFilename);
    if (!in.good()) {
        throw std::runtime_error("Nie można otworzyć pliku QUAL: " + qualFilename);
    }

    std::vector<std::vector<int>> allQuals;
    std::string line;
    std::vector<int> currentQuals;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // Nowa sekwencja
            if (!currentQuals.empty()) {
                allQuals.push_back(currentQuals);
            }
            currentQuals.clear();
        } else {
            // Czytamy wartości numeryczne jakości
            std::stringstream ss(line);
            int score;
            while (ss >> score) {
                currentQuals.push_back(score);
            }
        }
    }
    // Ostatnia sekwencja
    if (!currentQuals.empty()) {
        allQuals.push_back(currentQuals);
    }

    return allQuals;
}
