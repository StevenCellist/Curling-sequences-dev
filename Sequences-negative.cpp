// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 07-03-2021; expected time to length 48: 1 second* (yet to implement multi-threading)
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, Microsoft VS Studio 2019 Compiler

#include <iostream>
#include <vector>
#include <chrono>
#include <set>
#include <map>
#include <unordered_map>
#include <array>

#define PROFILING 
#ifdef PROFILING
#define PROFILE __declspec(noinline)
#else
#define PROFILE  
#endif
using namespace std::chrono;
const int length = 48;

int candidatecurl = 2;
int candidateperiod = 1;
std::vector<int> Tail = {};
std::vector<int> Periods = {};
std::vector<int> Generator = {};
std::vector<int> Max_tail_lengths = {};
std::map<int, std::vector<int>> Generators_memory = {};
std::vector<std::vector<int>> Best_generators = {};
std::set<int> Change_indices = {};

std::vector<int> seq_new = {};

std::unordered_map<int, int> expected_tails = {
    {2, 2},
    {4, 4},
    {6, 8},
    {8, 58},
    {9, 59},
    {10, 60},
    {11, 112},
    {14, 118},
    {19, 119},
    {22, 120},
    {48, 131},
    {68, 132},
    {73, 133},
    {77, 173}, // from http://neilsloane.com/doc/CNC.pdf
};

PROFILE
void krul(std::vector<int>* seq, int* curl, int* period) {           // curl = 1, period = 0
    int l = seq->size();
    for (int i = 1; i <= (l / 2); ++i) {
        int j = i;
        while ((*seq)[l - j - 1] == (*seq)[l - j - 1 + i]) {
            ++j;
            if (j >= l)
                break;
        }
        int candidate = j / i;
        if (candidate > *curl) {
            *curl = candidate;
            *period = i;
        }
    }
}

PROFILE
void tail_with_periods(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods) {               // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    krul(&seq, &curl, &period);
    while (curl > 1) {
        tail->push_back(curl);
        temp.push_back(curl);           // easy optimization
        periods->push_back(period);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
    }
}

PROFILE
void tail_with_periods_part(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods, int i) {   // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    krul(&seq, &curl, &period);
    while (curl > 1 and tail->size() < i) {
        tail->push_back(curl);
        temp.push_back(curl);           // easy optimization
        periods->push_back(period);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
    }
}

PROFILE
bool check_period_size() {
    return (candidatecurl * candidateperiod) > (length + Periods.size());
}

PROFILE
bool check_candidatecurl_size() {
    if (Tail.size())
        return candidatecurl > Tail.back() + 1;
    else
        return candidatecurl > length;
}

PROFILE
void up() {
    ++candidateperiod;
    while (check_period_size()) {
        ++candidatecurl;
        candidateperiod = 1;
        if (check_candidatecurl_size()) {
            auto index = std::find(Change_indices.begin(), Change_indices.end(), Periods.size());
            if (index != Change_indices.end()) {
                Change_indices.erase(index);
            }
            if (Tail.size() == 0) {
                candidatecurl = 0;
                break;
            }
            index = std::find(Change_indices.begin(), Change_indices.end(), Periods.size() - 1);
            if (index == Change_indices.end()) {
                Change_indices.insert(Periods.size() - 1);
                candidatecurl = Tail.back() + 1;
                candidateperiod = 1;
                Tail.pop_back();
                Periods.pop_back();
            }
            else {
                candidatecurl = Tail.back();
                candidateperiod = Periods.back() + 1;
                Tail.pop_back();
                Periods.pop_back();
                Generator = Generators_memory[Periods.size()];
                Generators_memory.erase(Periods.size());
            }
        }
    }
}

PROFILE
int real_generator_length() {
    int i = 0;
    while (Generator[i] == (-length + i)) {
        ++i;
        if (i == length)
            break;
    }
    return length - i;
}

PROFILE
bool check_positive(int len) {
    for (int i : std::vector<int>(Generator.begin() + length - len, Generator.end())) {
        if (i < 1)
            return false;
    }
    return true;
}

PROFILE
void append() {
    Generators_memory[Periods.size()] = Generator;
    Generator = std::vector<int>(seq_new.begin(), seq_new.begin() + length);

    Tail.push_back(candidatecurl);
    Periods.push_back(candidateperiod);

    int curl = 1, period = 0;
    std::vector<int> temp = Generator;
    temp.insert(temp.end(), Tail.begin(), Tail.end());
    while (true) {
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
        if (curl == 1)
            break;
        Tail.push_back(curl);
        temp.push_back(curl);
        Periods.push_back(period);
    }
    candidatecurl = 2;
    candidateperiod = 1;
    Change_indices.insert(Tail.size());
    int len = real_generator_length();
    if (Max_tail_lengths.back() == Tail.size()) {
        std::vector<int> temp = Best_generators.back();
        temp.insert(temp.end(), Generator.begin(), Generator.end() - len);
        Best_generators.back() = temp;
    }
    if (Max_tail_lengths[len - 1] < Tail.size()) {
        Max_tail_lengths[len - 1] = Tail.size();
        Best_generators[len - 1] = std::vector<int>(Generator.end() - len, Generator.end());
    }
}

PROFILE
bool test_1() {
    seq_new = Generator;
    seq_new.insert(seq_new.end(), Tail.begin(), Tail.end());

    int k = Generator.size();
    int l = seq_new.size();
    for (int i = 0; i < (candidatecurl - 1) * candidateperiod; ++i) {
        int a = seq_new[l - 1 - i];
        int b = seq_new[l - 1 - i - candidateperiod];
        if (a != b and a > 0 and b > 0)
            return false;
        if (a > b) {
            for (int j = 0; j < k; ++j) {
                if (seq_new[j] == b)
                    seq_new[j] = a;
            }
        }
        else if (a < b) {
            for (int j = 0; j < k; ++j) {
                if (seq_new[j] == a)
                    seq_new[j] = b;
            }
        }
    }
    return true;
}

PROFILE
bool test_2() {
    int l = seq_new.size();
    std::vector<int> temp_seq = seq_new;
    temp_seq.push_back(candidatecurl);
    std::vector<int> temp_period = Periods;
    temp_period.push_back(candidateperiod);
    int curl = 1, period = 0;
    for (int i = 0; i <= l - length; ++i) {
        std::vector<int> temp = std::vector<int>(seq_new.begin(), seq_new.begin() + length + i);
        curl = 1, period = 0;
        krul(&temp, &curl, &period);
        if (curl != temp_seq[length + i] or period != temp_period[i])
            return false;
    }
    return true;
}

PROFILE
bool check_if_period_works() {
    if (test_1()) {
        if (test_2())
            return true;
    }
    return false;
}

PROFILE
void backtracking_step() {
    if (check_if_period_works())
        append();
    else
        up();
}

PROFILE
void backtracking(int k1, int p1, int k2, int p2) {
    Change_indices.insert(0);
    for (int i = 0; i < length; ++i) {
        Generator.push_back(-length + i);
        Max_tail_lengths.push_back(0);
        Best_generators.push_back({});
    }
    std::vector<int> seq = Generator;
    seq.insert(seq.end(), Tail.begin(), Tail.end());

    while (candidatecurl) {
        if (Tail.size() == 0 and candidatecurl == k2 and candidateperiod == p2)
            break;
        backtracking_step();
    }
    int record = 0;
    for (int i = 0; i < length; ++i) {
        if (Max_tail_lengths[i] > record) {
            record = Max_tail_lengths[i];
            if (expected_tails.find(i + 1) == expected_tails.end())
                std::cout << "NEW:" << std::endl;
            else if (expected_tails[i + 1] != record)
                std::cout << "WRONG:" << std::endl;
            std::cout << i + 1 << ": " << record << ", [";
            for (int x : Best_generators[i])
                std::cout << x << ", ";
            std::cout << "]" << std::endl;
        }
    }
}

int main()
{
    //int k1, k2, p1, p2;
    //std::cout << "Length: ";
    //std::cin >> length;
    /*std::cin >> k1;
    std::cin >> p1;
    std::cin >> k2;
    std::cin >> p2;*/
    auto t1 = std::chrono::high_resolution_clock::now();
    backtracking(2, 1, 1000, 1000);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;

    //std::cout << "copies: " << copies << std::endl;
    //for (int i = 0; i < 250; ++i) {
    //    std::cout << i << "\trows:\t" << freq_row[i] << "\t, cols: \t" << freq_col[i] << std::endl;
    //}
}