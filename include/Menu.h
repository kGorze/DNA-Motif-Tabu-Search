#ifndef MENU_H
#define MENU_H

#include "EnhancedMotif.h"
#include <set>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "../include/FastaQualParser.h"
#include "../include/DNASequence.h"
#include "../include/MotifGraphBuilder.h"
#include "../include/TabuSearchMotif.h"
#include "../include/TabuSearchBase.h"
#include "../include/EnhancedMotif.h"

// Uruchamia główne menu programu
void runMenu();

// Nowa funkcja do uruchamiania automatycznych testów
void runAutomatedTests();

#endif // MENU_H
