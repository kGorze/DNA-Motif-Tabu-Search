cmake_minimum_required(VERSION 3.29)
project(zadanie_3)

set(CMAKE_CXX_STANDARD 14)

add_executable(zadanie_3
        source/main.cpp
        include/Graph.h
        source/Graph.cpp
        include/TabuSearchBase.h
        source/TabuSearchBase.cpp
        include/TabuSearchDeterministic.h
        source/TabuSearchDeterministic.cpp
        source/TabuSearchProbabilistic.cpp
        include/TabuSearchProbabilistic.h
        source/Stabulus.cpp
        include/Stabulus.h
        source/IteratedStabulus.cpp
        include/IteratedStabulus.h
        source/GraphGenerator.cpp
        include/GraphGenerator.h
        include/DNASequence.h
        source/DNASequence.cpp
        include/FastaQualParser.h
        source/FastaQualParser.cpp
        include/MotifGraphBuilder.h
        source/MotifGraphBuilder.cpp
        include/TabuSearchMotif.h
        source/TabuSearchMotif.cpp
        include/Menu.h
        source/Menu.cpp
        include/EnhancedMotif.h
        source/EnhancedMotif.cpp
)
