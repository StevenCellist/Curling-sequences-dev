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
#define FILE_OPEN file.open("Sequences.txt")
#define OUTPUT file
#define FILE_CLOSE file.close();
#endif

typedef std::vector<int16_t> val_vector;

const int length = 190;
const int thread_count = std::thread::hardware_concurrency();

struct context {
    int candidatecurl, candidateperiod;
    val_vector seq = val_vector(length), seq_new, periods, max_tails;
    std::vector<val_vector> best_generators;
    std::map<int16_t, val_vector> generators_memory;
    std::set<int16_t> change_indices = { 0 };
    std::array<std::vector<int>, 2 * length + 2> seq_map; // adjust for + 1 index on positive and negative side
    val_vector pairs;
    std::vector<int> temp;
};

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

int krul(const val_vector& s, int& period, int l, int minimum) {
    int curl = minimum - 1;                                 // base value for curl
    int limit = l / minimum;                                // limit up to which to check for repetition
    const int16_t* p1 = &s[l - 1];                      // start of the last pattern
    for (int i = 1; i <= limit; ++i, --p1) {                      // check for repetition up to the limit
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

__forceinline
void up(context& ctx) {
    ++ctx.candidateperiod;                                              // try period one larger now
    while (true) {
        if (ctx.periods.empty())                                        // stop if tail is empty
            break;
        if ((ctx.candidatecurl * ctx.candidateperiod) <= ctx.seq.size())        // if this pair is within sequence size, we can break and try that instead
            break;
        if (ctx.candidatecurl > ctx.seq.back()) {                           // if curl > last_curl, the sequence will have to change if it wants to become better
            int16_t k = (int16_t)ctx.periods.size() - 1;                // we will have to change the for-last curl in order to change the tail, 
            auto index = ctx.change_indices.find(k);                    // so let's check whether it has been done yet
            if (index == ctx.change_indices.end()) {                    // didn't find it?
                ctx.change_indices.insert(k);                           // insert it and let's try increased curl with new period
                ctx.candidatecurl = ctx.seq.back() + 1;
                ctx.candidateperiod = 1 + k / ctx.candidatecurl;
            }
            else {                                                  // did find it?
                ctx.candidatecurl = ctx.seq.back();                         // let's try this same curl but now with one higher period
                ctx.candidateperiod = ctx.periods.back() + 1;

                val_vector& temp = ctx.generators_memory[k];            // retrieve generator from memory
                for (int i = 0; i < length; i++) {
                    if (ctx.seq[i] != temp[i]) {                        // check where the differences are, apply them, and change map accordingly
                        for (int& x : ctx.seq_map[ctx.seq[i] + length])
                            if (x == i) {
                                x = ctx.seq_map[ctx.seq[i] + length].back();// use last value
                                break;
                            }
                        ctx.seq_map[ctx.seq[i] + length].pop_back();        // remove (now duplicated) last value
                        ctx.seq_map[temp[i] + length].push_back(i);
                        ctx.seq[i] = temp[i];
                    }
                }
                ctx.generators_memory.erase(k);                         // and delete memory entry
            }

            auto ind = std::find(ctx.seq_map[ctx.seq.back() + length].begin(), ctx.seq_map[ctx.seq.back() + length].end(), ctx.seq.size() - 1);
            ctx.seq_map[ctx.seq.back() + length].erase(ind);                // delete the last curl index from the map
            ctx.seq.pop_back();                                         // delete the last curl and its period from the tail, and delete entry from the changed indices
            ctx.periods.pop_back();

            index = ctx.change_indices.find(k + 1);
            if (index != ctx.change_indices.end())
                ctx.change_indices.erase(index);
        }
        else {
            ++ctx.candidatecurl;                                        // if not, we increase the curl to try and calculate its period
            ctx.candidateperiod = 1 + (int)ctx.periods.size() / ctx.candidatecurl;
        }
    }
}

int real_generator_length(context& ctx) {       // returns the index of the last unique value of the generator
    int i = 0;
    while ((ctx.seq[i] == (-length + i)) && (++i != length)) {}
    return (length - i);
}

__forceinline
void append(context& ctx) {
    for (int i = 0; i < ctx.pairs.size(); i += 2) {             // NOW we are sure we passed test_1 and test_2, so it is time to update the map
        for (int x : ctx.seq_map[ctx.pairs[i + 1] + length])
            ctx.seq_map[ctx.pairs[i] + length].push_back(x);        // so we move all changed values to their new location
        ctx.seq_map[ctx.pairs[i + 1] + length].clear();             // and delete their previous entries
    }
    ctx.seq.resize(length);                                     // discard the tail so we can swap it into the memory
    ctx.generators_memory[(int16_t)ctx.periods.size()].swap(ctx.seq);
    ctx.seq.swap(ctx.seq_new);                                      // retrieve the new sequence from test_2
    ctx.seq.push_back((int16_t)ctx.candidatecurl);                  // add the candidates because we passed test_1 and test_2
    ctx.periods.push_back((int16_t)ctx.candidateperiod);
    ctx.seq_map[ctx.candidatecurl + length].push_back(ctx.seq.size() - 1);
    int period = 0;
    while (true) {                                          // build the tail for this generator
        int curl = krul(ctx.seq, period, (int)ctx.seq.size(), 2);
        if (curl == 1)
            break;
        ctx.seq_map[curl + length].push_back(ctx.seq.size());
        ctx.seq.push_back((int16_t)curl);
        ctx.periods.push_back((int16_t)period);

    }

    int tail = (int)ctx.periods.size();                         // find tail size (seq.size() - length)
    ctx.candidatecurl = 2;                                      // prepare candidates for next backtracking_step
    ctx.candidateperiod = 1 + tail / 2;
    ctx.change_indices.insert((int16_t)tail);                   // to improve generator, we would have to take note of the index where the 1 occured
    int len = real_generator_length(ctx);                      // retrieve actual generator length
    if (ctx.max_tails[len] < (int16_t)tail) {                   // and update maximum values for this thread
        ctx.max_tails[len] = (int16_t)tail;
        ctx.best_generators[len] = val_vector(ctx.seq.begin() + length - len, ctx.seq.begin() + length);
    }
}

// this function checks whether the current sequence allows for the candidates to be added
__forceinline
bool test_1(context& ctx) {
    int l = (int)ctx.seq.size() - 1;                            // last element of pattern to check
    int lcp = l - ctx.candidateperiod;                          // last element of potential pattern
    int limit = (ctx.candidatecurl - 1) * ctx.candidateperiod;      // limit to which to check for repetition
    for (int i = 0; i < limit; ++i, --l, --lcp) {
        int16_t a = ctx.seq[l];
        int16_t b = ctx.seq[lcp];
        if (a != b and a > 0 and b > 0)                     // check whether the repetition may be possible
            return false;
    }
    ctx.seq_new = ctx.seq;                                          // create dummy sequence
    ctx.pairs.clear();
    l = (int)ctx.seq.size() - 1;                                // reset pattern values
    lcp = l - ctx.candidateperiod;

    for (int i = 0; i < limit; ++i, --l, --lcp) {
        int16_t a = ctx.seq_new[l];
        int16_t b = ctx.seq_new[lcp];
        if (a != b) {                                       // because we are changing values below, we may encounter a possible break, again
            if (a > 0 and b > 0)
                return false;
            if (a > b)
                std::swap(a, b);                            // a is now always < b and < 0
            ctx.pairs.push_back(b);                             // add (b, a) combo to the pairs
            ctx.pairs.push_back(a);

            ctx.temp.clear();
            ctx.temp.push_back(a);                                  // temporary vector that will hold all map values that need to be changed
            for (int index = 0; index < ctx.temp.size(); index++)   // because we don't (want to) change the map here (not sure we pass test_1 and test_2)
                for (int i = 0; i < ctx.pairs.size(); i += 2)       // if we change a to b, and later change b, we also need to change a in that case
                    if (ctx.pairs[i] == ctx.temp[index])                // so we need to check if we already crossed the value b
                        ctx.temp.push_back(ctx.pairs[i + 1]);
            for (int x : ctx.temp)                              // apply changes to sequence
                for (int ind : ctx.seq_map[x + length])
                    ctx.seq_new[ind] = b;
        }
    }
    return true;
}

// this function checks whether the proposed change invalidates the generator (regarding curl or period)
__forceinline
bool test_2(context& ctx) {
    int l = (int)ctx.seq_new.size();
    int period = 0;

    for (int i = 0; i < l - length; ++i) {                                  // check within tail for a valid curl or period
        int curl = krul(ctx.seq_new, period, length + i, ctx.seq_new[length + i]);  // calculate curl and period up to this part of the sequence
        if (curl != ctx.seq_new[length + i] or period != ctx.periods[i])            // if the curl or its period are not related, the tail will not improve
            return false;
    }
    int curl = krul(ctx.seq_new, period, l, ctx.candidatecurl);
    return (curl == ctx.candidatecurl and period == ctx.candidateperiod);
}

void backtracking_step(context& ctx) {
    if (test_1(ctx) && test_2(ctx))   // depending on whether the sequence will improve the tail length...
        append(ctx);               // we make the tail longer if it passed the checks
    else
        up(ctx);                   // or we upgrade the generator if they failed
}

void backtracking() {
    context ctx;
    for (int i = 0; i <= length; ++i) {                  // initiate thread-local record vectors
        ctx.max_tails.push_back(0);
        ctx.best_generators.push_back({});
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    while (true) {
        {
            std::lock_guard<std::mutex> l(m_kp);        // lock g_k1 and g_p1
            if (g_k1 > length)
                break;
            ctx.candidatecurl = g_k1;                       // value for first curl of the tail
            ctx.candidateperiod = g_p1;                     // value for period of this curl
            int limit = length / g_k1;
            if (++g_p1 > limit) {                       // gone out of range?
                g_p1 = 1;                               // reset period
                ++g_k1;                                 // try new curl
            }
        }
        for (int i = 0; i < length; ++i)                // initiate sequence
            ctx.seq[i] = (int16_t)(-length + i);

        //seq_map.clear();
        //seq_map.resize(2 * length + 2);
        for (auto& v : ctx.seq_map)
            v.clear();
        for (int j = 0; j < length; ++j)
            ctx.seq_map[ctx.seq[j] + length].push_back(j);

        backtracking_step(ctx);                            // perform backtracking for this combination (g_k1, g_p1)
        while (ctx.periods.size())
            backtracking_step(ctx);
    }
    {
        std::lock_guard<std::mutex> l(m_tails);         // update global arrays under a safe lock
        for (int i = 0; i <= length; ++i) {
            if (ctx.max_tails[i] > g_max_tails[i]) {
                g_max_tails[i] = ctx.max_tails[i];
                g_best_generators[i] = ctx.best_generators[i];
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
}