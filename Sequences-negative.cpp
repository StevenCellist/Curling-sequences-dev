// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 07-03-2021; expected time to length 48: 1 second*
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

#ifdef _MSC_VER
#define FILE_OPEN   
#define OUTPUT std::cout
#define FILE_CLOSE  
#ifdef PROFILING
#define PROFILE __declspec(noinline)
#else
#define PROFILE  
#endif
#else // gcc (Linux)
#include <fstream>
std::ofstream file;
#define FILE_OPEN file.open("Sequences.txt")
#define OUTPUT file
#define FILE_CLOSE file.close();
#define PROFILE  
#endif

typedef int16_t val_type;
typedef std::vector<val_type> val_vector;
using namespace std::chrono;

const int length = 80;

thread_local int candidatecurl, candidateperiod;
thread_local val_vector seq(length), seq_new, Periods, Max_tails;
thread_local std::vector<val_vector> Best_generators;
thread_local std::map<val_type, val_vector> Generators_memory = {};
thread_local std::set<val_type> Change_indices = { 0 };

val_vector Global_max_tails(length, 0);
std::vector<val_vector> Global_best_generators(length);
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
    {149, 343},
    {154, 356},
    {176, 406},
    {197, 1668},
    {198, 1669},
    {200, 1670}
};

bool diff(const val_type* p1, const val_type* p2, int count) {
    while (count--) {
        if (*p1++ != *p2++)
            return true;
    }
    return false;
}

PROFILE
int krul(const val_vector& s, int& period, int l) {
    int curl = 1;                                           // base value for curl
    period = 0;                                             // this is redundant but does speed up by 4%(?!)
    int limit = l / 2;                                      // limit up to which to check for repetition
    for (int i = 1; i <= limit; ++i) {                      // check for repetition up to the limit
        const val_type* p1 = &s[l - i];                     // start of the last pattern
        const val_type* p2 = p1 - i;                        // start of the previous pattern
        for (int freq = 2; p2 >= &s[0]; ++freq, p2 -= i) {
            if (diff(p1, p2, i))                            // doesn't match
                break;
            if (curl < freq) {                              // found better curl?
                curl = freq;                                // update curl
                limit = l / (curl + 1);                     // update limit
                period = i;                                 // update period
            }
        }
    }
    return curl;
}

bool check_candidatecurl_size() {
    if (seq.size() > length)
        return candidatecurl > seq.back() + 1;
    else
        return candidatecurl > length;
}

PROFILE
void up() {
    ++candidateperiod;
    while ((candidatecurl * candidateperiod) > seq.size()) {
        ++candidatecurl;
        candidateperiod = 1 + Periods.size() / candidatecurl;
        if (check_candidatecurl_size()) {
            if (seq.size() == length) {
                candidatecurl = 0;
                break;
            }
            val_type k = (val_type)Periods.size();
            auto index = Change_indices.find(k);
            if (index != Change_indices.end())
                Change_indices.erase(index);
            index = Change_indices.find(k - 1);
            if (index == Change_indices.end()) {
                Change_indices.insert(k - 1);
                candidatecurl = seq.back() + 1;
                candidateperiod = 1 + (k - 1) / candidatecurl;
            }
            else {
                candidatecurl = seq.back();
                candidateperiod = Periods.back() + 1;
                memcpy(&seq[0], &Generators_memory[k - 1][0], length * sizeof(val_type));
                Generators_memory.erase(k - 1);
            }
            seq.pop_back();
            Periods.pop_back();
        }
    }
}

int real_generator_length() {       // returns the index of the last unique value of the generator
    int i = 0;
    while ((seq[i] == (-length + i)) && (++i != length)) {}
    return (length - i);
}

PROFILE
void append() {
    seq.resize(length);
    Generators_memory[(val_type)Periods.size()].swap(seq);

    seq.swap(seq_new);
    Periods.push_back((val_type)candidateperiod);

    int period = 0;
    while (true) {
        int curl = krul(seq, period, (int)seq.size());
        if (curl == 1)
            break;
        seq.push_back((val_type)curl);
        Periods.push_back((val_type)period);
    }

    candidatecurl = 2;
    int tail = (int)Periods.size();
    candidateperiod = 1 + tail / 2;
    Change_indices.insert((val_type)tail);
    int len = real_generator_length();
    if (Max_tails[len - 1] < (val_type)tail) {
        Max_tails[len - 1] = (val_type)tail;
        Best_generators[len - 1] = val_vector(seq.begin() + length - len, seq.begin() + length);
    }
}

PROFILE
bool test_1() {

    int l = (int)seq.size() - 1;
    int lcp = l - candidateperiod;
    int limit = (candidatecurl - 1) * candidateperiod;
    for (int i = 0; i < limit; ++i, --l, --lcp) {
        val_type a = seq[l];
        val_type b = seq[lcp];
        if (a != b and a > 0 and b > 0)
            return false;
    }
    seq_new = seq;
    l = (int)seq.size() - 1;
    lcp = l - candidateperiod;

    int jmax = length - 1;      // get the index of the last negative value in seq_new to limit the for loop below
    for (; jmax > 0; --jmax) {
        if (seq_new[jmax] < 0)
            break;
    }

    for (int i = 0; i < limit; ++i, --l, --lcp) {
        val_type a = seq_new[l];
        val_type b = seq_new[lcp];
        if (a != b) {
            if (a > 0 and b > 0)
                return false;
            if (a > b)
                std::swap(a, b); // a is now always < b
            for (int j = 0; j <= jmax; ++j) {   // don't need to loop through positive values.
                if (seq_new[j] == a) {
                    seq_new[j] = b;
                    if (j == jmax)
                        for (--jmax; jmax > 0; --jmax) {
                            if (seq_new[jmax] < 0)
                                break;
                        }
                }
            }
        }
    }
    return true;
}

PROFILE
bool test_2() {
    int l = (int)seq_new.size();
    seq_new.push_back((val_type)candidatecurl);
    val_vector temp_period = Periods;
    temp_period.push_back((val_type)candidateperiod);
    int period = 0;

    for (int i = 0; i <= l - length; ++i) {
        int curl = krul(seq_new, period, length + i);
        if (curl != seq_new[length + i] or period != temp_period[i])
            return false;
    }
    return true;
}

PROFILE
void backtracking_step() {
    if (test_1() && test_2())   // depending on whether the sequence will improve the tail length...
        append();               // we make the tail longer if it passed the checks
    else
        up();                   // or we upgrade the generator if they failed
}

PROFILE
void backtracking(std::vector<int> params) {
    candidatecurl = params[0];                          // start-value for first curl of the tail
    candidateperiod = params[1];                        // start-value for period of first curl of the tail
    int k2 = params[2];                                 // stop-value for first curl of the tail
    int p2 = params[3];                                 // stop-value for period of first curl of the tail

    for (int i = 0; i < length; ++i) {                  // initiate sequence and its thread-local companions
        seq[i] = (val_type)(-length + i);
        Max_tails.push_back(0);
        Best_generators.push_back({});
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    while (candidatecurl) {                             // perform backtracking while we are within thread value boundaries
        if (seq.size() == length and candidatecurl == k2 and candidateperiod == p2)
            break;
        backtracking_step();
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    {
        std::lock_guard<std::mutex> l(m_tails);         // update global arrays under a safe lock
        for (int i = 0; i < length; ++i) {
            if (Max_tails[i] > Global_max_tails[i]) {
                Global_max_tails[i] = Max_tails[i];
                Global_best_generators[i] = Best_generators[i];
            }
        }
        OUTPUT << "Finished: " << length << ", " << params[0] << ", " << params[1] << ", " << k2 << ", " << p2 << ", ";
        OUTPUT << "duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;
    }
}

void multi_threader() {
    std::vector<std::thread> thread_vector;
    std::vector<std::vector<int>> params;   // parameters per thread
    if constexpr (length <= 80)             // for easy debugging and fast benchmarking
        params = { {2, 1, 1000, 1000} };
    else if constexpr (length < 110)        // for small multi-threaded runs
        params = { {2, 1, 2, 3}, {2, 3, 2, 7}, {2, 7, 2, 24}, {2, 24, 2, 40}, {2, 40, 3, 3}, {3, 3, 3, 24}, {3, 24, 5, 1}, {5, 1, 1000, 1000} };
    else                                    // 8- or 16-threaded parameters for long runs
        params = { {2, 1, 2, 3}, {2, 3, 2, 6}, {2, 6, 2, 20}, {2, 20, 2, 60}, {2, 60, 3, 2}, {3, 2, 4, 1}, {4, 1, 5, 1}, {5, 1, 1000, 1000} };
    //params = { {2, 1, 2, 2}, {2, 2, 2, 4}, {2, 4, 2, 7}, {2, 7, 2, 12}, {2, 12, 2, 25}, {2, 25, 2, 60}, {2, 60, 2, 100}, {2, 100, 3, 1},
    //           {3, 1, 3, 3}, {3, 3, 3, 20}, {3, 20, 3, 60}, {3, 60, 4, 2}, {4, 2, 5, 1}, {5, 1, 6, 3}, {6, 2, 7, 1}, {7, 1, 1000, 1000}};
    for (std::vector<int> x : params) 
        thread_vector.emplace_back(std::thread(backtracking, x));
    for (auto& th : thread_vector) 
        th.join();
}

int main()
{
    FILE_OPEN;
    auto t1 = std::chrono::high_resolution_clock::now();
    multi_threader();
    auto t2 = std::chrono::high_resolution_clock::now();

    int record = 0;
    for (int i = 0; i < length; ++i) {
        if (Global_max_tails[i] > record) {
            record = Global_max_tails[i];
            if (expected_tails.find(i + 1) == expected_tails.end())
                OUTPUT << "NEW:" << std::endl;
            else if (expected_tails[i + 1] != record)
                OUTPUT << "WRONG:" << std::endl;
            OUTPUT << i + 1 << ": " << record << ", [";
            for (int x : Global_best_generators[i])
                OUTPUT << x << ",";
            OUTPUT << "]" << std::endl;
        }
    }
    FILE_CLOSE;
}