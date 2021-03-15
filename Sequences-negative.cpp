// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 15-03-2021; completed time to length 48: 0.5 second*
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
#include <cstring>
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
const int thread_count = std::thread::hardware_concurrency();

thread_local int candidatecurl, candidateperiod;
thread_local val_vector seq(length), seq_new, Periods, Max_tails;
thread_local std::vector<val_vector> Best_generators;
thread_local std::map<val_type, val_vector> Generators_memory = {};
thread_local std::set<val_type> Change_indices = { 0 };

val_vector Global_max_tails(length + 1, 0);
std::vector<val_vector> Global_best_generators(length + 1);
std::mutex m_tails, m_kp;

int k1 = 2, p1 = 1;
const int k2 = 1000, p2 = 1000;

std::unordered_map<int, int> expected_tails = {
    {2, 2},     {4, 4},     {6, 8},     {8, 58},    {9, 59},     {10, 60},    {11, 112},   {14, 118},   {19, 119},  {22, 120}, 
    {48, 131},  {68, 132},  {73, 133},  {77, 173},  {85, 179},   {115, 215},  {116, 228},  {118, 229},  {128, 332}, {132, 340},
    {133, 342}, {149, 343}, {154, 356}, {176, 406}, {197, 1668}, {198, 1669}, {200, 1670}, {208, 1706}, {217, 1836}
};

bool diff(const val_type* p1, const val_type* p2, int count) {
    while (count--) {
        if (*p1++ != *p2++)
            return true;
    }
    return false;
}

PROFILE
int krul(const val_vector& s, int& period, int l, int minimum) {
    int curl = minimum - 1;                                 // base value for curl
    int limit = l / minimum;                                // limit up to which to check for repetition
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

PROFILE
void up() {
    ++candidateperiod;                                              // try period one larger now
    while (true) {
        if (!Periods.size())                                        // stop if tail is empty
            break;                                 
        if ((candidatecurl * candidateperiod) <= seq.size())        // if this pair is within sequence size, we can break and try that instead
            break; 
        if (candidatecurl > seq.back()) {                           // if curl > last_curl, the sequence will have to change if it wants to become better
            val_type k = (val_type)Periods.size() - 1;              // we will have to change the for-last curl in order to change the tail, 
            auto index = Change_indices.find(k);                    // so let's check whether it has been done yet
            if (index == Change_indices.end()) {                    // didn't find it?
                Change_indices.insert(k);                           // insert it and let's try increased curl with new period
                candidatecurl = seq.back() + 1;
                candidateperiod = 1 + k / candidatecurl;
            }
            else {                                                  // did find it?
                candidatecurl = seq.back();                         // let's try this same curl but now with one higher period
                candidateperiod = Periods.back() + 1;
                memcpy(&seq[0], &Generators_memory[k][0], length * sizeof(val_type)); // retrieve the original generator from memory
                Generators_memory.erase(k);                         // and delete it
            }
            seq.pop_back();                                         // delete the last curl and its period from the tail, and delete entry from the changed indices
            Periods.pop_back();
            index = Change_indices.find(k + 1);
            if (index != Change_indices.end())
                Change_indices.erase(index);
        }
        else {
            ++candidatecurl;                                        // if not, we increase the curl to try and calculate its period
            candidateperiod = 1 + (int)Periods.size() / candidatecurl;
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
    seq.resize(length);                                     // discard the tail so we can swap it into the memory
    Generators_memory[(val_type)Periods.size()].swap(seq);
    seq.swap(seq_new);                                      // retrieve the new sequence from test_2
    Periods.push_back((val_type)candidateperiod);

    int period = 0;
    while (true) {                                          // build the tail for this generator
        int curl = krul(seq, period, (int)seq.size(), 2);
        if (curl == 1)
            break;
        seq.push_back((val_type)curl);
        Periods.push_back((val_type)period);
    }

    int tail = (int)Periods.size();                         // find tail size (seq.size() - length)
    candidatecurl = 2;                                      // prepare candidates for next backtracking_step
    candidateperiod = 1 + tail / 2;
    Change_indices.insert((val_type)tail);                  // to improve generator, we would have to take note of the index where the 1 occured
    int len = real_generator_length();                      // retrieve actual generator length
    if (Max_tails[len] < (val_type)tail) {                  // and update maximum values for this thread
        Max_tails[len] = (val_type)tail;
        Best_generators[len] = val_vector(seq.begin() + length - len, seq.begin() + length);
    }
}

PROFILE
// this function checks whether the current sequence allows for the candidates to be added
bool test_1() {
    int l = (int)seq.size() - 1;                            // last element of pattern to check
    int lcp = l - candidateperiod;                          // last element of potential pattern
    int limit = (candidatecurl - 1) * candidateperiod;      // limit to which to check for repetition
    for (int i = 0; i < limit; ++i, --l, --lcp) {
        val_type a = seq[l];
        val_type b = seq[lcp];
        if (a != b and a > 0 and b > 0)                     // check whether the repetition may be possible
            return false;
    }
    seq_new = seq;                                          // create dummy sequence
    l = (int)seq.size() - 1;                                // reset pattern values
    lcp = l - candidateperiod;

    int jmax = length - 1;                                  // get the index of the last negative value in seq_new to limit the for loop below
    for (; jmax > 0; --jmax) {
        if (seq_new[jmax] < 0)
            break;
    }
    for (int i = 0; i < limit; ++i, --l, --lcp) {
        val_type a = seq_new[l];
        val_type b = seq_new[lcp];
        if (a != b) {                                       // because we are changing values below, we may encounter a possible break, again
            if (a > 0 and b > 0)
                return false;
            if (a > b)
                std::swap(a, b);                            // a is now always < b and < 0
            for (int j = 0; j <= jmax; ++j) {               // don't need to loop through positive values.
                if (seq_new[j] == a) {
                    seq_new[j] = b;
                    if (j == jmax)                          // check whether we can break earlier
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
// this function checks whether the tail will contain anything meaningful
bool test_2() {
    int l = (int)seq_new.size();
    seq_new.push_back((val_type)candidatecurl);                             // add the candidates because we passed test_1
    int period = 0;

    for (int i = 0; i < l - length; ++i) {                                  // check within tail for a valid curl or period
        int curl = krul(seq_new, period, length + i, seq_new[length + i]);  // calculate curl and period up to this part of the sequence
        if (curl != seq_new[length + i] or period != Periods[i])            // if the curl or its period are not related, the tail will not improve
            return false;
    }
    int curl = krul(seq_new, period, l, seq_new[l]);
    return (curl == seq_new[l] and period == candidateperiod);
}

PROFILE
void backtracking_step() {
    if (test_1() && test_2())   // depending on whether the sequence will improve the tail length...
        append();               // we make the tail longer if it passed the checks
    else
        up();                   // or we upgrade the generator if they failed
}

PROFILE
void backtracking() {
    for (int i = 0; i <= length; ++i) {                  // initiate thread-local record vectors
        Max_tails.push_back(0);
        Best_generators.push_back({});
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    while (true) {
        {
            std::lock_guard<std::mutex> l(m_kp);        // lock k1 and p1
            if (k1 == k2)
                break;
            candidatecurl = k1;                         // value for first curl of the tail
            candidateperiod = p1;                       // value for period of this curl
            if (++p1 == p2) {
                p1 = 1;
                ++k1;
            }
        }
        for (int i = 0; i < length; ++i)                // initiate sequence
            seq[i] = (val_type)(-length + i);

        backtracking_step();                            // perform backtracking for this combination (k1, p1)
        while (Periods.size())
            backtracking_step();
    }
    {
        std::lock_guard<std::mutex> l(m_tails);         // update global arrays under a safe lock
        for (int i = 0; i <= length; ++i) {
            if (Max_tails[i] > Global_max_tails[i]) {
                Global_max_tails[i] = Max_tails[i];
                Global_best_generators[i] = Best_generators[i];
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    OUTPUT << "Finished: " << length << ", " << "duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;
}

void multi_threader() {
    std::vector<std::thread> thread_vector;
    if (length <= 80)
        thread_vector.emplace_back(std::thread(backtracking));
    else 
        for (int i = 0; i < thread_count; ++i)
            thread_vector.emplace_back(std::thread(backtracking));
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
    for (int i = 0; i <= length; ++i) {
        if (Global_max_tails[i] > record) {
            record = Global_max_tails[i];
            if (expected_tails.find(i) == expected_tails.end())
                OUTPUT << "NEW:" << std::endl;
            else if (expected_tails[i] != record)
                OUTPUT << "WRONG:" << std::endl;
            OUTPUT << i << ": " << record << ", [";
            for (int x : Global_best_generators[i])
                OUTPUT << x << ",";
            OUTPUT << "]" << std::endl;
        }
    }
    FILE_CLOSE;
}