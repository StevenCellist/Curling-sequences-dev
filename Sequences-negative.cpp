// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 19-03-2021; completed time to length 48: 98 milliseconds*
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

#ifdef _MSC_VER
#define FILE_OPEN   
#define OUTPUT std::cout
#define FILE_CLOSE  
#else // gcc (Linux)
#include <fstream>
#include <cstring>
std::ofstream file;
#define FILE_OPEN file.open("Sequences.txt", ios_base::out | ios_base::app)
#define OUTPUT file
#define FILE_CLOSE file.close();
#endif

typedef std::vector<int16_t> val_vector;

const int length = 160;
const int thread_count = std::thread::hardware_concurrency();

thread_local int candidatecurl, candidateperiod;
thread_local val_vector seq(length), seq_new, periods, max_tails;
thread_local std::vector<val_vector> best_generators;
thread_local std::map<int16_t, val_vector> generators_memory = {};
thread_local std::set<int16_t> change_indices = { 0 };
thread_local uint64_t c1(0), c2(0);
uint64_t g_c1(0), g_c2(0);

val_vector g_max_tails(length + 1, 0);
std::vector<val_vector> g_best_generators(length + 1);
std::mutex m_tails, m_kp;

int g_k1 = 2, g_p1 = 1;

std::unordered_map<int, int> expected_tails = {
    {2, 2},     {4, 4},     {6, 8},     {8, 58},    {9, 59},     {10, 60},    {11, 112},   {14, 118},   {19, 119},   {22, 120},
    {48, 131},  {68, 132},  {73, 133},  {77, 173},  {85, 179},   {115, 215},  {116, 228},  {118, 229},  {128, 332},  {132, 340},
    {133, 342}, {149, 343}, {154, 356}, {176, 406}, {197, 1668}, {199, 1669}, {200, 1670}, {208, 1708}, {217, 1836}, {290, 3382}
};

bool diff(const int16_t* p1, const int16_t* p2, int count) {
    int count64 = count / 4;    // four 16-bit ints in uint64_t
    if (count64) {
        uint64_t* p1_64 = (uint64_t*)p1;
        uint64_t* p2_64 = (uint64_t*)p2;
        while (count64--) {
            if (*p1_64++ != *p2_64++)
                return true;
        }
    }

    int count16 = count % 4;    // left-overs
    if (count16) {
        int offset = count - count16;
        p1 += offset;
        p2 += offset;
        while (count16--) {
            if (*p1++ != *p2++)
                return true;
        }
    }
    return false;
}

int krul(const val_vector& s, int& period, int l) {
    int curl = 1;                                           // base value for curl
    int limit = l / 2;                                      // limit up to which to check for repetition
    for (int i = 1; i <= limit; ++i) {                      // check for repetition up to the limit
        const int16_t* p1 = &s[l - i];                      // start of the last pattern
        const int16_t* p2 = p1 - i;                         // start of the previous pattern
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

void up() {
    c1++;
    ++candidateperiod;                                              // try period one larger now
    while (true) {
        if (periods.empty())                                        // stop if tail is empty
            break;
        if ((candidatecurl * candidateperiod) <= seq.size())        // if this pair is within sequence size, we can break and try that instead
            break;
        if (candidatecurl > seq.back()) {                           // if curl > last_curl, the sequence will have to change if it wants to become better
            int16_t k = (int16_t)periods.size() - 1;                // we will have to change the for-last curl in order to change the tail, 
            auto index = change_indices.find(k);                    // so let's check whether it has been done yet
            if (index == change_indices.end()) {                    // didn't find it?
                change_indices.insert(k);                           // insert it and let's try increased curl with new period
                candidatecurl = seq.back() + 1;
                candidateperiod = 1 + k / candidatecurl;
            }
            else {                                                  // did find it?
                candidatecurl = seq.back();                         // let's try this same curl but now with one higher period
                candidateperiod = periods.back() + 1;
                memcpy(&seq[0], &generators_memory[k][0], length * sizeof(int16_t)); // retrieve the original generator from memory
                c2++;
                generators_memory.erase(k);                         // and delete it
            }
            seq.pop_back();                                         // delete the last curl and its period from the tail, and delete entry from the changed indices
            periods.pop_back();
            index = change_indices.find(k + 1);
            if (index != change_indices.end())
                change_indices.erase(index);
        }
        else {
            ++candidatecurl;                                        // if not, we increase the curl to try and calculate its period
            candidateperiod = 1 + (int)periods.size() / candidatecurl;
        }
    }
}

int real_generator_length() {       // returns the index of the last unique value of the generator
    int i = 0;
    while ((seq[i] == (-length + i)) && (++i != length)) {}
    return (length - i);
}

void append() {
    seq.resize(length);                                     // discard the tail so we can swap it into the memory
    generators_memory[(int16_t)periods.size()].swap(seq);
    seq.swap(seq_new);                                      // retrieve the new sequence from test_2
    periods.push_back((int16_t)candidateperiod);

    int period = 0;
    while (true) {                                          // build the tail for this generator
        int curl = krul(seq, period, (int)seq.size());
        if (curl == 1)
            break;
        seq.push_back((int16_t)curl);
        periods.push_back((int16_t)period);
    }

    int tail = (int)periods.size();                         // find tail size (seq.size() - length)
    candidatecurl = 2;                                      // prepare candidates for next backtracking_step
    candidateperiod = 1 + tail / 2;
    change_indices.insert((int16_t)tail);                   // to improve generator, we would have to take note of the index where the 1 occured
    int len = real_generator_length();                      // retrieve actual generator length
    if (max_tails[len] < (int16_t)tail) {                   // and update maximum values for this thread
        max_tails[len] = (int16_t)tail;
        best_generators[len] = val_vector(seq.begin() + length - len, seq.begin() + length);
    }
}

int get_first_index(const val_vector& s) {                  // returns the index of the first modofied value of the generator
    int i = 0;
    while ((s[i] == (-length + i)) && (++i != length)) {}
    return i;
}

// this function checks whether the current sequence allows for the candidates to be added
bool test_1() {
    int l = (int)seq.size() - 1;                            // last element of pattern to check
    int lcp = l - candidateperiod;                          // last element of potential pattern
    int limit = (candidatecurl - 1) * candidateperiod;      // limit to which to check for repetition
    for (int i = 0; i < limit; ++i, --l, --lcp) {
        int16_t a = seq[l];
        int16_t b = seq[lcp];
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
    int first = get_first_index(seq_new);

    for (int i = 0; i < limit; ++i, --l, --lcp) {
        int16_t a = seq_new[l];
        int16_t b = seq_new[lcp];
        if (a != b) {                                       // because we are changing values below, we may encounter a possible break, again
            if (a > 0 and b > 0)
                return false;
            if (a > b)
                std::swap(a, b);                            // a is now always < b and < 0
            for (int j = std::min(first, lcp); j <= jmax; ++j) {    // don't need to loop through positive values at the end or unique at the start
                if (seq_new[j] == a) {
                    seq_new[j] = b;
                    if (j == jmax)                          // check whether we can break earlier
                        for (--jmax; jmax > first; --jmax) {
                            if (seq_new[jmax] < 0)
                                break;
                        }
                }
            }
        }
    }
    return true;
}

// this function checks whether the proposed change invalidates the generator (regarding curl or period)
bool test_2() {
    int l = (int)seq_new.size();
    seq_new.push_back((int16_t)candidatecurl);                              // add the candidates because we passed test_1
    int period = 0;

    for (int i = 0; i < l - length; ++i) {                                  // check within tail for a valid curl or period
        int curl = krul(seq_new, period, length + i);                       // calculate curl and period up to this part of the sequence
        if (curl != seq_new[length + i] or period != periods[i])            // if the curl or its period are not related, the tail will not improve
            return false;
    }
    int curl = krul(seq_new, period, l);
    return (curl == seq_new[l] and period == candidateperiod);
}

void backtracking_step() {
    if (test_1() && test_2())   // depending on whether the sequence will improve the tail length...
        append();               // we make the tail longer if it passed the checks
    else
        up();                   // or we upgrade the generator if they failed
}

void backtracking() {
    for (int i = 0; i <= length; ++i) {                  // initiate thread-local record vectors
        max_tails.push_back(0);
        best_generators.push_back({});
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    while (true) {
        {
            std::lock_guard<std::mutex> l(m_kp);        // lock g_k1 and g_p1
            if (g_k1 > length)
                break;
            candidatecurl = g_k1;                       // value for first curl of the tail
            candidateperiod = g_p1;                     // value for period of this curl
            int limit = length / g_k1;
            if (++g_p1 > limit) {                       // gone out of range?
                g_p1 = 1;                               // reset period
                ++g_k1;                                 // try new curl
            }
        }
        for (int i = 0; i < length; ++i)                // initiate sequence
            seq[i] = (int16_t)(-length + i);

        backtracking_step();                            // perform backtracking for this combination (g_k1, g_p1)
        while (periods.size())
            backtracking_step();
    }
    {
        std::lock_guard<std::mutex> l(m_tails);         // update global arrays under a safe lock
        for (int i = 0; i <= length; ++i) {
            if (max_tails[i] > g_max_tails[i]) {
                g_max_tails[i] = max_tails[i];
                g_best_generators[i] = best_generators[i];
            }
        }
        g_c1 += c1;
        g_c2 += c2;
        std::cout << "c1: " << c1 << ", c2: " << c2 << std::endl;
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
        if (g_max_tails[i] > record) {
            record = g_max_tails[i];
            if (expected_tails.find(i) == expected_tails.end())
                OUTPUT << "NEW:" << std::endl;
            else if (expected_tails[i] != record)
                OUTPUT << "WRONG:" << std::endl;
            OUTPUT << i << ": " << record << ", [";
            for (int x : g_best_generators[i])
                OUTPUT << x << ",";
            OUTPUT << "]" << std::endl;
        }
    }
    FILE_CLOSE;
    std::cout << "c1: " << g_c1 << ", c2: " << g_c2 << std::endl;
}