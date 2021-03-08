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
#include <thread>
#include <mutex>

#define PROFILING 
#ifdef PROFILING
#define PROFILE __declspec(noinline)
#else
#define PROFILE  
#endif
using namespace std::chrono;
int length = 48;

thread_local int candidatecurl;
thread_local int candidateperiod;
thread_local std::vector<int> Tail = {};
thread_local std::vector<int> Periods = {};
thread_local std::vector<int> Generator = {};
thread_local std::vector<int> Max_tails = {};
thread_local std::vector<int> seq_new = {};
thread_local std::vector<std::vector<int>> Best_generators = {};
thread_local std::map<int, std::vector<int>> Generators_memory = {};
thread_local std::set<int> Change_indices = {};

std::vector<int> Global_max_tails(256, 0);
std::vector<std::vector<int>> Global_best_generators(256);
std::mutex m_tails;

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
    {77, 173},
    {85, 179},
    {115, 215},
    {116, 228},
    {118, 229},
    {128, 332},
    {132, 340},
    {133, 342},
    {149, 343}
};

PROFILE
void krul(std::vector<int>* seq, int* curl, int* period) {           // curl = 1, period = 0
    *curl = 1, * period = 0;
    int l = seq->size();
    for (int i = 1; i <= (l / 2); ++i) {
        int j = i;
        while ((*seq)[l - j - 1] == (*seq)[l - j - 1 + i]) {
            ++j;
            if (j >= l) {
                int candidate = j / i;
                if (candidate > *curl) {
                    *curl = candidate;
                    *period = i;
                }
                break;
            }
            int candidate = j / i;
            if (candidate > *curl) {
                *curl = candidate;
                *period = i;
            }
        }
    }
}

PROFILE
void tail_with_periods(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods) {               // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    int l = seq.size();
    krul(&seq, &curl, &period);
    while (curl > 1) {
        temp.push_back(curl);
        periods->push_back(period);
        krul(&temp, &curl, &period);
    }
    *tail = std::vector<int>(temp.begin() + l, temp.end());
}

PROFILE
void tail_with_periods_part(std::vector<int> seq, std::vector<int>* tail, std::vector<int>* periods, int i) {   // tail = {}, periods = {}
    int curl = 1, period = 0;
    std::vector<int> temp = seq;
    int l = seq.size();
    krul(&seq, &curl, &period);
    while (curl > 1 and tail->size() < i) {
        temp.push_back(curl);
        periods->push_back(period);
        krul(&temp, &curl, &period);
    }
    *tail = std::vector<int>(temp.begin() + l, temp.end());
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
            if (!Tail.size()) {
                candidatecurl = 0;
                break;
            }
            int k = Periods.size();
            auto index = std::find(Change_indices.begin(), Change_indices.end(), k);
            if (index != Change_indices.end()) {
                Change_indices.erase(index);
            }
            index = std::find(Change_indices.begin(), Change_indices.end(), k - 1);
            if (index == Change_indices.end()) {
                Change_indices.insert(k - 1);
                candidatecurl = Tail.back() + 1;
                candidateperiod = 1;
            }
            else {
                candidatecurl = Tail.back();
                candidateperiod = Periods.back() + 1;
                Generator = Generators_memory[k - 1];
                Generators_memory.erase(k - 1);
            }
            Tail.pop_back();
            Periods.pop_back();
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
    if (Max_tails.back() == Tail.size()) {
        std::vector<int> temp = Best_generators.back();
        temp.insert(temp.end(), Generator.begin(), Generator.end() - len);
        Best_generators.back() = temp;
    }
    if (Max_tails[len - 1] < Tail.size()) {
        Max_tails[len - 1] = Tail.size();
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
    candidatecurl = k1;
    candidateperiod = p1;
    Change_indices.insert(0);
    for (int i = 0; i < length; ++i) {
        Generator.push_back(-length + i);
        Max_tails.push_back(0);
        Best_generators.push_back({});
    }
    std::vector<int> seq = Generator;
    seq.insert(seq.end(), Tail.begin(), Tail.end());

    auto t1 = std::chrono::high_resolution_clock::now();
    while (candidatecurl) {
        if (Tail.size() == 0 and candidatecurl == k2 and candidateperiod == p2)
            break;
        backtracking_step();
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    {
        std::lock_guard<std::mutex> l(m_tails);
        for (int i = 0; i < length; ++i) {
            if (Max_tails[i] > Global_max_tails[i]) {
                Global_max_tails[i] = Max_tails[i];
                Global_best_generators[i] = Best_generators[i];
            }
        }
        std::cout << "Finished: " << k1 << ", " << p1 << ", " << k2 << ", " << p2 << ", ";
        std::cout << "duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;
    }
}

void multi_threader() {
    std::vector<std::thread> thread_vector;
    // good values for length > 110
    /*thread_vector.emplace_back(std::thread(backtracking, 2, 1, 2, 3));
    thread_vector.emplace_back(std::thread(backtracking, 2, 3, 2, 6));
    thread_vector.emplace_back(std::thread(backtracking, 2, 6, 2, 20));
    thread_vector.emplace_back(std::thread(backtracking, 2, 20, 2, 60));
    thread_vector.emplace_back(std::thread(backtracking, 2, 60, 3, 2));
    thread_vector.emplace_back(std::thread(backtracking, 3, 2, 4, 1));
    thread_vector.emplace_back(std::thread(backtracking, 4, 1, 5, 1));
    thread_vector.emplace_back(std::thread(backtracking, 5, 1, 1000, 1000));
    // good values for length 80 ~ 110
    thread_vector.emplace_back(std::thread(backtracking, 2, 1, 2, 3));
    thread_vector.emplace_back(std::thread(backtracking, 2, 3, 2, 7));
    thread_vector.emplace_back(std::thread(backtracking, 2, 7, 2, 24));
    thread_vector.emplace_back(std::thread(backtracking, 2, 24, 2, 40));
    thread_vector.emplace_back(std::thread(backtracking, 2, 40, 3, 3));
    thread_vector.emplace_back(std::thread(backtracking, 3, 3, 3, 24));
    thread_vector.emplace_back(std::thread(backtracking, 3, 24, 5, 1));
    thread_vector.emplace_back(std::thread(backtracking, 5, 1, 1000, 1000));*/
    // no point in multi-threading for length < 80
    thread_vector.emplace_back(std::thread(backtracking, 2, 1, 1000, 1000));
    for (auto& th : thread_vector) th.join();
}

int main()
{
    /*for (int i = 50; i <= 120; i += 5) {
        length = i;
        auto t1 = std::chrono::high_resolution_clock::now();
        multi_threader();
        auto t2 = std::chrono::high_resolution_clock::now();
        int record = 0;
        for (int i = 0; i < length; ++i) {
            if (Global_max_tails[i] > record) {
                record = Global_max_tails[i];
                if (expected_tails.find(i + 1) == expected_tails.end())
                    std::cout << "NEW:" << std::endl;
                else if (expected_tails[i + 1] != record)
                    std::cout << "WRONG:" << std::endl;
                std::cout << i + 1 << ": " << record << ", [";
                for (int x : Global_best_generators[i])
                    std::cout << x << ", ";
                std::cout << "]" << std::endl;
            }
        }
        std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;
    }*/
    auto t1 = std::chrono::high_resolution_clock::now();
    multi_threader();
    auto t2 = std::chrono::high_resolution_clock::now();

    int record = 0;
    for (int i = 0; i < length; ++i) {
        if (Global_max_tails[i] > record) {
            record = Global_max_tails[i];
            if (expected_tails.find(i + 1) == expected_tails.end())
                std::cout << "NEW:" << std::endl;
            else if (expected_tails[i + 1] != record)
                std::cout << "WRONG:" << std::endl;
            std::cout << i + 1 << ": " << record << ", [";
            for (int x : Global_best_generators[i])
                std::cout << x << ", ";
            std::cout << "]" << std::endl;
        }
    }
    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;

    //std::cout << "copies: " << copies << std::endl;
    //for (int i = 0; i < 250; ++i) {
    //    std::cout << i << "\trows:\t" << freq_row[i] << "\t, cols: \t" << freq_col[i] << std::endl;
    //}
}