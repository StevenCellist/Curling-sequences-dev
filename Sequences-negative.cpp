// Sequences-negative.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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
const int length = 68;

//std::array<uint64_t, 250> freq_row = {}, freq_col = {};
//uint64_t copies = 0;

class TDict {
public:
    PROFILE
        TDict() {
        //std::memset(&m_table[0][0], -1, m_rows * m_cols * sizeof int16_t);
    }
    PROFILE
        void insert(int ind, int16_t value) {
        int row = ind + m_bias;
        if (m_row_max < row)
            m_row_max = row;

        int col_max = m_col_max[row];
        for (int col = 0; col < col_max; ++col) {
            if (m_table[row][col] == value) {
                return; // was already there, nothing to insert
            }
        }
        m_table[row][col_max] = value;
        m_col_max[row]++; // one value added
    }
    PROFILE
        void erase(int ind, int16_t value) {
        int row = ind + m_bias;
        int col_max = m_col_max[row];
        for (int col = 0; col <= col_max; ++col) {
            if (m_table[row][col] == value) { // found the value to remove
                m_col_max[ind + m_bias]--; // one value removed
                for (; col < col_max; ++col) { // shift the rest of the row left by 1
                    m_table[row][col] = m_table[row][col + 1];
                }
                break;
            }
        }
    }
    PROFILE
        void erase(int ind) {
        m_col_max[ind + m_bias] = 0; // all values removed
    }
    PROFILE
        void move(int src, int dst) {
        if (m_row_max < dst + m_bias)
            m_row_max = dst + m_bias;
        int col_max = m_col_max[src + m_bias];
        for (int col = 0; col < col_max; ++col) {
            insert(dst, m_table[src + m_bias][col]);
        }

        erase(src);
    }
    PROFILE
        bool is_size_1(int src) {
        return m_col_max[src + m_bias] == 1;
        //return m_table[src + m_bias][0] != -1 && m_table[src + m_bias][1] == -1;
    }

    // copy assignment operator; could use the default, but this is a bit faster as it only copies up to that.m_row_max
    PROFILE
        TDict& operator=(const TDict& that) {
        //copies++;
        if (this != &that) {
            //freq_row[that.m_row_max]++;
            for (int row = 0; row <= that.m_row_max; ++row) {
                int col_max = that.m_col_max[row];
                //freq_col[col_max]++;
                //if(col_max == 1)
                //    m_table[row][0] = that.m_table[row][0];
                //else if (col_max > 1)
                    for (int col = 0; col < col_max; ++col) {
                        m_table[row][col] = that.m_table[row][col];
                    }
            }

            int row_max = std::max(m_row_max, that.m_row_max);
            std::memcpy(&m_col_max, &that.m_col_max, (row_max + 1) * sizeof int16_t);
            m_row_max = that.m_row_max;
            //m_col_max = that.m_col_max;
        }
        return *this;
    }
private:
    static const int m_bias = length-1; // adjust to 0-based index; we index from [-(length-1)]
    static const int m_rows = 250;      // TODO: Steven - is this a `length + max tail`?
    static const int m_cols = 250;      // TODO: Vlad - rework to `vectors`, so that we don't have to fix the array length
    std::array<std::array<int16_t, m_cols>, m_rows> m_table;
    int m_row_max = 0;
    std::array<int16_t, m_rows> m_col_max = {};
};

int candidatecurl = 2;
int candidateperiod = 1;
std::vector<int> Tail = {};
std::vector<int> Periods = {};
std::vector<int> Generator = {};
std::vector<int> Max_tail_lengths = {};
std::map<int, std::vector<int>> Generators_memory = {};
std::vector<std::vector<int>> Best_generators = {};
std::set<int> Change_indices = {};
TDict Dict;
TDict Dict_new;
std::map<int, TDict> Dicts_memory = {};

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
    {76, 133},
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
                //Dict[Tail.back()].erase(length + Tail.size() - 1);
                Dict.erase(Tail.back(), length + Tail.size() - 1);
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
                Dict = Dicts_memory[Periods.size()];
                Generators_memory.erase(Periods.size());
                Dicts_memory.erase(Periods.size());
            }
        }
    }
}

PROFILE
int real_generator_length() {
    int i = 0;
    while (Dict.is_size_1(Generator[i])) {
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
    Dicts_memory[Periods.size()] = Dict;
    Generator = std::vector<int>(seq_new.begin(), seq_new.begin() + length);
    Dict = Dict_new;
    Dict.insert(candidatecurl, length + Tail.size());
    Dict_new.insert(candidatecurl, length + Tail.size());

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
        Dict.insert(curl, length + Tail.size() - 1);
        Dict_new.insert(curl, length + Tail.size() - 1);
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

    Dict_new = Dict;

    int l = seq_new.size();
    for (int i = 0; i < (candidatecurl - 1) * candidateperiod; ++i) {
        int a = seq_new[l - 1 - i];
        int b = seq_new[l - 1 - i - candidateperiod];
        if (a != b and a > 0 and b > 0)
            return false;
        if (a > b) {
            for (int j = 0; j < l; ++j) {
                if (seq_new[j] == b)
                    seq_new[j] = a;
            }
            Dict_new.move(b, a);
        }
        else if (a < b) {
            for (int j = 0; j < l; ++j) {
                if (seq_new[j] == a)
                    seq_new[j] = b;
            }
            Dict_new.move(a, b);
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
        Generator.push_back(-i);
        Dict.insert(-i, i);
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
            if(expected_tails.find(i+1) == expected_tails.end())
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
