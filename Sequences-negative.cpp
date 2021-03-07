// Sequences-negative.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <set>
#include <map>
#include <unordered_map>
#include <array>

using namespace std::chrono;
const int length = 68;

class TDict {
public:
    __declspec(noinline)
        TDict() {
        std::memset(&m_table[0][0], -1, m_rows * m_cols * sizeof int16_t);
    }
    __declspec(noinline)
        void insert(int ind, int16_t value) {
        if (m_row_max < ind + m_bias)
            m_row_max = ind + m_bias;

        for (auto& x : m_table[ind + m_bias]) {
            if (x == value) {
                break;
            }
            if (x == -1) {
                x = value;
                break;
            }
        }
    }
    __declspec(noinline)
        void erase(int ind, int16_t value) {
        int i = 0;
        for (; i < m_table[ind + m_bias].size(); ++i) {
            if (m_table[ind + m_bias][i] == -1) {
                return; // was not there, nothing to erase
            }
            if (m_table[ind + m_bias][i] == value) {
                break;
            }
        }
        for (; i < m_table[ind + m_bias].size() - 1; ++i) {
            if ((m_table[ind + m_bias][i] = m_table[ind + m_bias][i + 1]) == -1) {
                return; // just copied the sentinel
            }
        }
    }
    __declspec(noinline)
        void erase(int ind) {
        std::memset(&m_table[ind + m_bias], -1, m_cols * sizeof int16_t);
    }
    __declspec(noinline)
        void move(int src, int dst) {
        if (m_row_max < dst + m_bias)
            m_row_max = dst + m_bias;
        for (auto& x : m_table[src + m_bias]) {
            if (x == -1)
                break;
            insert(dst, x);
        }
        erase(src);
    }
    __declspec(noinline)
        bool is_size_1(int src) {
        return m_table[src + m_bias][0] != -1 && m_table[src + m_bias][1] == -1;
    }

    // copy assignment operator; could use the default, but this is a bit faster as it only copies up to that.m_row_max
    __declspec(noinline)
        TDict& operator=(const TDict& that) {
        if (this != &that) {
            for (int row = 0; row <= that.m_row_max; ++row) {
                for (int col = 0; col < m_cols; ++col) {
                    if ((m_table[row][col] = that.m_table[row][col]) == -1) {
                        break; // just copied the sentinel
                    }
                }
            }

            for (int row = that.m_row_max + 1; row <= m_row_max; ++row) {// clear not copied rows
                std::memset(&m_table[row], -1, m_cols * sizeof int16_t);
            }
            m_row_max = that.m_row_max;

        }
        return *this;
    }
private:
    static const int m_bias = length-1; // adjust to 0-based index; we index from [-(length-1)]
    static const int m_rows = 250;      // TODO: Steven - is this a `length + max tail`?
    static const int m_cols = 250;      // TODO: Vlad - rework to `vectors`, so that we don't have to fix the array length
    std::array<std::array<int16_t, m_cols>, m_rows> m_table;
    int m_row_max = 0;
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

__declspec(noinline)
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

__declspec(noinline)
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

__declspec(noinline)
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

__declspec(noinline)
bool check_period_size() {
    return (candidatecurl * candidateperiod) > (length + Periods.size());
}

__declspec(noinline)
bool check_candidatecurl_size() {
    if (Tail.size())
        return candidatecurl > Tail.back() + 1;
    else
        return candidatecurl > length;
}

__declspec(noinline)
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

__declspec(noinline)
int real_generator_length() {
    int i = 0;
    while (Dict.is_size_1(Generator[i])) {
        ++i;
        if (i == length)
            break;
    }
    return length - i;
}

__declspec(noinline)
bool check_positive(int len) {
    for (int i : std::vector<int>(Generator.begin() + length - len, Generator.end())) {
        if (i < 1)
            return false;
    }
    return true;
}

__declspec(noinline)
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

__declspec(noinline)
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

__declspec(noinline)
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

__declspec(noinline)
bool check_if_period_works() {
    if (test_1()) {
        if (test_2())
            return true;
    }
    return false;
}

__declspec(noinline)
void backtracking_step() {
    if (check_if_period_works())
        append();
    else
        up();
}

__declspec(noinline)
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
}
