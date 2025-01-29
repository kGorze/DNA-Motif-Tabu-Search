# Wykrywanie Motywów DNA za pomocą Tabu Search

## Spis treści
- [Wprowadzenie](#wprowadzenie)
- [Struktura projektu](#struktura-projektu)
- [Algorytmy i metody](#algorytmy-i-metody)
- [Główne komponenty](#główne-komponenty)
- [Konfiguracja](#konfiguracja)
- [Interfejsy i klasy bazowe](#interfejsy-i-klasy-bazowe)
- [Szczegółowa dokumentacja](#szczegółowa-dokumentacja)
- [Format plików wejściowych/wyjściowych](#format-plików-wejściowychwyjściowych)
- [Sposób użycia](#sposób-użycia)
- [Wymagania systemowe](#wymagania-systemowe)
- [Testy i weryfikacja](#testy-i-weryfikacja)

## Wprowadzenie
Niniejszy projekt zawiera implementację zaawansowanych technik wyszukiwania i analizy **motywów w sekwencjach DNA** przy wykorzystaniu podejścia **Tabu Search** (TS). Motyw DNA to krótki fragment (tzw. k-mer) pojawiający się w kilku sekwencjach DNA z pewną liczbą dopuszczalnych błędów (mismatchy). Celem algorytmu jest odnalezienie wspólnego fragmentu (lub grupy spójnych fragmentów) dla wybranego zbioru sekwencji.

W projekcie zaimplementowano zarówno **deterministyczne** (np. *TabuSearchDeterministic*) jak i **probabilistyczne** (np. *TabuSearchProbabilistic*) warianty metody TS, a także wyspecjalizowaną wersję do wykrywania motywów w sekwencjach biologicznych (*TabuSearchMotif*, *EnhancedTabuSearchMotif*).

## Struktura projektu

```
.
├── include/                           
│   ├── DNASequence.h                 # Reprezentacja pojedynczej sekwencji DNA
│   ├── EnhancedMotif.h               # Rozszerzone struktury i klasy do motywów
│   ├── FastaQualParser.h             # Parser plików FASTA/QUAL
│   ├── Graph.h                       # Struktura grafu i dane o wierzchołkach
│   ├── GraphGenerator.h              # Generator losowych grafów
│   ├── IteratedStabulus.h            # Iterowany algorytm do znajdowania zbiorów stabilnych
│   ├── Menu.h                        # Menu interaktywne do obsługi programu
│   ├── MotifGraphBuilder.h           # Klasa do budowania grafu motywów
│   ├── Stabulus.h                    # Prosta klasa do znajdowania zbioru stabilnego
│   ├── TabuSearchBase.h              # Klasa bazowa Tabu Search
│   ├── TabuSearchDeterministic.h     # Deterministyczna implementacja TS
│   ├── TabuSearchMotif.h             # Specjalizacja TS dla problemu motywów
│   ├── TabuSearchProbabilistic.h     # Probabilistyczna implementacja TS
│   └── ...
├── source/
│   ├── DNASequence.cpp               # Implementacja klasy DNASequence
│   ├── EnhancedMotif.cpp             # Implementacja rozszerzonych klas motywów
│   ├── FastaQualParser.cpp           # Implementacja parsera FASTA/QUAL
│   ├── Graph.cpp                     # Implementacja struktur grafowych
│   ├── GraphGenerator.cpp            # Implementacja generatora grafów
│   ├── IteratedStabulus.cpp         # Implementacja algorytmu iterowanego
│   ├── Menu.cpp                      # Implementacja menu i logiki interakcji
│   ├── MotifGraphBuilder.cpp         # Logika budowania grafu dla detekcji motywów
│   ├── Stabulus.cpp                  # Algorytm do zbiorów stabilnych
│   ├── TabuSearchBase.cpp            # Metody wspólne dla wszystkich wariantów TS
│   ├── TabuSearchDeterministic.cpp   # Deterministyczny TS
│   ├── TabuSearchMotif.cpp           # Tabu Search wyspecjalizowany dla motywów
│   ├── TabuSearchProbabilistic.cpp   # Probabilistyczny TS
│   └── main.cpp                      # Główny punkt wejścia programu
└── CMakeLists.txt                    # Plik konfiguracyjny CMake
```

## Algorytmy i metody

### 1. Tabu Search (TS)
Podstawowy schemat wyszukiwania:
1. **Inicjalizacja** – tworzenie rozwiązania startowego (np. heurystyka zachłanna).
2. **Definicja sąsiedztwa** – zbiór możliwych przekształceń (dodanie/usuniecie wierzchołka z kliki).
3. **Listy tabu** – przechowywanie informacji o zakazanych ruchach, ograniczenie zapętlania.
4. **Kryterium aspiracji** – zezwalanie na ruch tabu, jeśli poprawia najlepsze dotąd rozwiązanie.
5. **Kryterium zatrzymania** – liczba iteracji bez poprawy, limit czasu, itp.

### 2. Deterministyczny vs. Probabilistyczny TS
- **TabuSearchDeterministic** – za każdym razem wybierany jest *najlepszy* kandydat z sąsiedztwa (zgodnie z ustalonym rankingiem).
- **TabuSearchProbabilistic** – próbkowanie sąsiedztwa, losowe wstrząsy (shake-up) i wybór najlepszego kandydata z próby.

### 3. Wyszukiwanie motywów w DNA
- **TabuSearchMotif** – wariant TS dostosowany do problemu motywów, z dodatkowymi funkcjami:
  - *evaluateMotifSolution(...)* – ocenia "jakość" motywu, biorąc pod uwagę wielkość kliki, liczbę sekwencji, różnice pomiędzy k-merami i ich lokalizację w sekwencjach.
  - *validMotifSolution(...)* – weryfikuje, czy rozwiązanie odpowiada poprawnemu motywowi (np. każdy wierzchołek w klastrze odpowiada innej sekwencji DNA).

- **EnhancedTabuSearchMotif** (w `EnhancedMotif.cpp/.h`) – rozszerza powyższy algorytm o dodatkowe kryteria (np. entropia, "average Hamming distance" itp.) oraz mechanizmy dywersyfikacji.

## Główne komponenty

### 1. `DNASequence`
Reprezentuje jedną sekwencję DNA wraz z informacjami o jakości (jeśli są dostępne):
```cpp
class DNASequence {
public:
    std::string name;                // Nazwa sekwencji (nagłówek FASTA)
    std::string bases;               // Ciąg nukleotydów (A, C, G, T)
    std::vector<int> qualityScores;  // Ścieżka jakości (opcjonalnie)
};
```

### 2. `FastaQualParser`
Wczytuje pary plików `.fasta` i `.qual`, tworząc obiekty `DNASequence`. 
- Metoda `parseFastaAndQual(fastaFilename, qualFilename)` zwraca wektor sekwencji.

### 3. `Graph` i `VertexData`
W projekcie motyw interpretuje się jako klikę w grafie:
```cpp
struct VertexData {
    int sequenceIndex;  // Indeks sekwencji (w wektorze DNASequence)
    int position;       // Pozycja startowa k-mera w sekwencji
    std::string kmer;   // Sam tekst k-mera
    int quality;        // (opcjonalnie) Średnia jakość
};

class Graph {
public:
    int n;
    std::vector<std::vector<bool>> adj;
    std::vector<VertexData> vertexInfo;

    Graph(int n_);
    void add_edge(int u, int v);
    bool is_edge(int u, int v) const;
    ...
};
```

**Koncepcja**: Wierzchołek = potencjalny k-mer, krawędź łączy dwa k-mery, jeśli spełniają kryterium podobieństwa i pozycja jest "zgodna" z parametrami (np. *allowedMismatches*, *positionMultiplier*).

### 4. `MotifGraphBuilder` / `EnhancedMotifGraphBuilder`
- Budują graf motywów z podanego zbioru sekwencji.
- Filtrowanie k-merów na podstawie jakości (`qualityThreshold`), długości k-min / k-max itp.

### 5. `Stabulus` / `IteratedStabulus`
- Proste algorytmy do wyszukiwania stabilnych zbiorów (i.e. niezależnych) w grafie, używane w niektórych eksperymentach (lub jako narzędzia pomocnicze).

### 6. `TabuSearchBase` i specjalizacje
Zawiera wspólne elementy Tabu Search:
```cpp
struct Solution {
    std::vector<bool> inClique;  // Czy dany wierzchołek należy do rozwiązania
    int size;                    // Rozmiar kliki
};

class TabuSearchBase {
protected:
    const Graph &G;
    ...
    virtual void initialize() = 0;
    virtual Solution selectBestNeighbor(const Solution &S,
                                        const std::vector<Solution> &neighbors) = 0;
    ...
public:
    virtual Solution run() = 0;
};
```

**Specjalizacje**:
- `TabuSearchDeterministic`
- `TabuSearchProbabilistic`
- `TabuSearchMotif`
- `EnhancedTabuSearchMotif` (z **Enhanced** elementami jak entropia motywu, adaptacja itp.)

## Konfiguracja
Projekt używa **CMake**. Główne kroki kompilacji:

1. Utwórz katalog do budowania:
   ```bash
   mkdir build
   cd build
   ```
2. Uruchom CMake:
   ```bash
   cmake ..
   ```
3. Skompiluj:
   ```bash
   make
   ```

Parametry algorytmu (np. `kMin, kMax, allowedMismatches, qualityThreshold`) wprowadzane są zwykle interaktywnie przez menu lub przez dopisanie w kodzie w `main.cpp` / `Menu.cpp`.

## Interfejsy i klasy bazowe

### `TabuSearchBase`
Główne metody:
- `run()` – punkt wejścia do algorytmu.
- `initialize()` – inicjalizacja (heurystyka startowa).
- `selectBestNeighbor(...)` – wybór najlepszego sąsiada według wybranej strategii.
- `updateTabuListAfterMove(...)` – aktualizacja listy tabu (T1, T2).
  
### `TabuSearchDeterministic` / `TabuSearchProbabilistic`
Dziedziczą po `TabuSearchBase` i implementują specyficzne strategie:
- **Deterministyczna**: szuka *najlepszego* ruchu w całym sąsiedztwie.
- **Probabilistyczna**: próbkowanie i/lub losowe perturbacje (shake-up), by uniknąć lokalnych maksimów.

### `TabuSearchMotif`
Dziedziczy po `TabuSearchBase`, zawiera:
- `validMotifSolution(...)` – sprawdza, czy rozwiązanie to poprawny motyw (jedna k-mera na sekwencję).
- `evaluateMotifSolution(...)` – oblicza ocenę na bazie podobieństwa k-merów i innych kryteriów.

### `EnhancedTabuSearchMotif`
Rozszerza `TabuSearchMotif` o:
- Metryki jakości (entropia, Hamming distance)
- Zaawansowaną dywersyfikację (np. adaptacyjne usuwanie wierzchołków)
- Dodatkowe logowanie i raporty.

## Szczegółowa dokumentacja

### Plik `main.cpp` i `Menu.cpp`
- `runMenu()` – interaktywne menu:
  1. Wczytaj sekwencje DNA z plików `sample.fasta` / `sample.qual` (przykład).
  2. Wybór 5 sekwencji, wstrzyknięcie motywu (*inject motif*) z określonymi parametrami (modyfikacje w sekwencjach).
  3. Budowa grafu k-merów (`MotifGraphBuilder` lub `EnhancedMotifGraphBuilder`).
  4. Uruchomienie Tabu Search w trybie *motif* (dobór parametrów: `kMin`, `kMax`, `qualityThreshold`, `allowedMismatches` itp.).
  5. Eksperymenty automatyczne (np. `runAutomatedTests()`), generowanie plików CSV z wynikami.

### Klasa `EnhancedMotifGraphBuilder`
- Rozszerza budowanie grafu o dodatkową weryfikację wstrzykniętych motywów, debug informację itp.

### Klasa `EnhancedTabuSearchMotif`
- Nadpisuje `evaluateMotifSolution(...)` tak, by uwzględniać zróżnicowane kryteria (m.in. entropia, różnorodność, penalizacje za rozbieżność pozycji k-merów).
- Zawiera logikę dywersyfikacji w metodzie `diversify(...)` dostosowaną do długotrwałych stagnacji.

## Format plików wejściowych/wyjściowych

### 1. Pliki FASTA/QUAL
- **FASTA** (`.fasta`): sekwencje w formacie:
  ```
  >Seq1
  ACGTACGTACGT...
  >Seq2
  GTTACGTCGATA...
  ```
- **QUAL** (`.qual`): analogicznie, pliki z oceną jakości:
  ```
  >Seq1
  40 39 39 40 38 ...
  >Seq2
  35 35 40 40 37 ...
  ```
Parser `FastaQualParser` kojarzy linie w pliku `.qual` z odpowiadającymi im sekwencjami z pliku `.fasta`.

### 2. Plik `LOAD_INSTANCE.txt`
- Przykład formatu stanu instancji (zapisywana/odczytywana przez `saveCurrentTestInstanceToFile()` i `loadTestInstanceFromFile()` w `Menu.cpp`):
  ```
  INSTANCE 1
  MOTIF ACGT
  SEQUENCES 5
  SEQ_START
  NAME MySeq1
  LENGTH 100
  BASES AAA...ACGT...TTT
  QUALITIES 39 40 [40 40 40] 40 ...
  MOTIF_POS 25
  SEQ_END
  ...
  ```

### 3. Raporty CSV
- Generowane przez `runAutomatedTests()`, zawierają informacje o czasie wykonania, wielkości znalezionych motywów, jakości, entropii, itp.

## Sposób użycia

1. **Kompilacja**:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```
2. **Uruchomienie**:
   ```bash
   ./dna_motif_tabu
   ```
3. **Interaktywne menu**:
   - Wczytaj pliki FASTA/QUAL.
   - Wybierz 5 sekwencji (lub załaduj predefiniowane `testX.fasta`/`testX.qual`).
   - Dodaj (zainjektuj) motyw o zadanej długości 4..9.
   - Zbuduj graf i uruchom TabuSearchMotif z ustalonymi parametrami.

4. **Tryb automatycznych testów**:
   - Opcja `runAutomatedTests()` generuje raport CSV i testuje różne progi jakości, różne długości motywu itp.

## Wymagania systemowe

- **Kompilator C++17** (lub nowszy)
- **CMake 3.10** (lub nowszy)
- (Opcjonalnie) **Boost / GoogleTest** – do testów jednostkowych
- System Linux/Windows/macOS (sprawdzany głównie na Linux)

## Testy i weryfikacja

- **Testy jednostkowe** (GoogleTest / Catch2) mogą obejmować:
  - Poprawność wczytywania plików FASTA/QUAL.
  - Działanie `TabuSearchBase` i jego metod (np. zarządzanie listami tabu).
  - Weryfikację budowy grafu z `MotifGraphBuilder`.
- **Testy integracyjne** sprawdzają współpracę algorytmu TS z modułem budowy grafu i parserem sekwencji.
- **Testy wydajnościowe** (stress test) – np. uruchamianie programu na dużych sekwencjach i pomiar czasu.
- **Raporty** w plikach CSV – weryfikacja, czy parametry i wyniki są zgodne z oczekiwaniami.

Przykład skryptu uruchamiającego testy (Linux):
```bash
#!/bin/bash
cd build
make tests
./tests/unit_tests
```

## Autorzy
- Konrad Gorzelańczyk  

## Data
styczeń 2025
