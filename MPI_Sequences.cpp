// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 11-06-2021; estimated time to length 48: 68 milliseconds*
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, Microsoft VS Studio 2019 Compiler

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <algorithm>
#include <mpi.h>

#define FILE_OPEN file.open("Sequences.txt")
#define OUTPUT file
#define FILE_CLOSE file.close();
#define INLINING inline __attribute__((always_inline))

typedef std::vector<int16_t> v16_t;
std::ofstream file;

const int length = 100;     // Tweakable parameter: set this to the desired generator length (n)
const int g_limit2 = 16;    // Tweakable parameter: increase this value if ranks do not finish simultaneously (necessary for large # of ranks)

struct context {                                            // all necessary variables for a rank
    bool c2p2 = false, c3p3 = false;
    int c_cand = 0, p_cand = 0, max_tails[length + 1] = { 0 };
    v16_t seq = v16_t(length), seq_new, periods, pairs, temp;
    int16_t best_generators[length + 1][length] = { 0 };
    std::map<int16_t, v16_t> generators_memory;
    std::unordered_set<int16_t> change_indices = { 0 };
    std::array<std::vector<int>, 2 * length + 2> seq_map;   // adjust for + 1 index on positive and negative side
};

std::unordered_map<int, int> expected_tails = {             // the record tail lengths we found so far
    {2, 2},     {4, 4},     {6, 8},     {8, 58},    {9, 59},     {10, 60},    {11, 112},   {14, 118},   {19, 119},   {22, 120},
    {48, 131},  {68, 132},  {73, 133},  {77, 173},  {85, 179},   {115, 215},  {116, 228},  {118, 229},  {128, 332},  {132, 340},
    {133, 342}, {149, 343}, {154, 356}, {176, 406}, {197, 1668}, {199, 1669}, {200, 1670}, {208, 1708}, {217, 1836}, {290, 3382}
};

bool diff(const int16_t* p1, const int16_t* p2, int count) {
    int count64 = count / 4;    // compare four 16-bit ints in uint64_t simultaneously
    if (count64) {
        uint64_t* p1_64 = (uint64_t*)p1;
        uint64_t* p2_64 = (uint64_t*)p2;
        while (count64--) {
            if (*p1_64++ != *p2_64++)
                return true;
        }
        p1 = (int16_t*)p1_64;
        p2 = (int16_t*)p2_64;
    }

    int count16 = count % 4;    // left-overs
    if (count16) {
        while (count16--) {
            if (*p1++ != *p2++)
                return true;
        }
    }
    return false;
}

// find the curl of the sequence
int krul(const v16_t& s, int& period, int l, int minimum) {
    int curl = minimum - 1;                                 // base value for curl
    int limit = l / minimum;                                // limit up to which to check for repetition
    const int16_t* p1 = &s[l - 1];                          // start of the last pattern
    for (int i = 1; i <= limit; ++i, --p1) {                // check for repetition up to the limit
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

// erase element from vector
void erase(std::vector<int>& v, int x) {
    int i = 0;
    while (v[i] != x)
        ++i;
    v[i] = v.back();        // use last value
    v.pop_back();           // remove (now duplicated) last value
}

INLINING // modify the candidate(s) and maybe the sequence to try and maximize the tail length
void up(context& ctx) {
    ++ctx.p_cand;                                                           // try period one larger now
    while (true) {
        if (ctx.periods.size() == ctx.c2p2 + ctx.c3p3)                      // stop if tail is empty (relative to depth)
            break;
        if ((ctx.c_cand * ctx.p_cand) <= ctx.seq.size())                    // if this pair is within sequence size, we can break and try that instead
            break;
        if (ctx.c_cand > ctx.seq.back()) {                                  // if curl > last_curl, the sequence will have to change to become longer
            int16_t k = (int16_t)ctx.periods.size() - 1;                    // we will have to change the for-last curl in order to change the tail, 
            auto index = ctx.change_indices.find(k);                        // so let's check whether it has been done yet
            if (index == ctx.change_indices.end()) {                        // didn't find it?
                ctx.change_indices.insert(k);                               // insert it and let's try increased curl with new period
                ctx.c_cand = ctx.seq.back() + 1;
                ctx.p_cand = 1 + k / ctx.c_cand;
            }
            else {                                                          // did find it?
                ctx.c_cand = ctx.seq.back();                                // let's try this same curl but now with one higher period
                ctx.p_cand = ctx.periods.back() + 1;

                v16_t& temp = ctx.generators_memory[k];                     // retrieve generator from memory
                for (int i = 0; i < length; i++) {
                    if (ctx.seq[i] != temp[i]) {                            // check where the differences are, apply them, and change map accordingly
                        erase(ctx.seq_map[ctx.seq[i] + length], i);
                        ctx.seq_map[temp[i] + length].push_back(i);
                        ctx.seq[i] = temp[i];
                    }
                }
                ctx.generators_memory.erase(k);                             // and delete memory entry
            }

            auto ind = std::find(ctx.seq_map[ctx.seq.back() + length].begin(), ctx.seq_map[ctx.seq.back() + length].end(), ctx.seq.size() - 1);
            ctx.seq_map[ctx.seq.back() + length].erase(ind);                // delete the last curl, period and index from their vectors
            ctx.seq.pop_back();
            ctx.periods.pop_back();

            index = ctx.change_indices.find(k + 1);
            if (index != ctx.change_indices.end())
                ctx.change_indices.erase(index);
        }
        else {
            ++ctx.c_cand;                                                   // if not, we increase the curl to try and calculate its period
            ctx.p_cand = 1 + (int)ctx.periods.size() / ctx.c_cand;
        }
    }
}

// find the shortest generator that could have generated this tail
int real_generator_length(context& ctx) {
    int i = 0;
    while ((ctx.seq[i] == (-length + i)) && (++i != length)) {}
    return (length - i);
}

INLINING // construct the tail for this generator
void append(context& ctx) {
    for (int i = 0; i < ctx.pairs.size(); i += 2) {             // NOW we are sure we passed test_1 and test_2, so it is time to update the map
        for (int x : ctx.seq_map[ctx.pairs[i + 1] + length])
            ctx.seq_map[ctx.pairs[i] + length].push_back(x);    // so we move all changed values to their new location
        ctx.seq_map[ctx.pairs[i + 1] + length].clear();         // and delete their previous entries
    }
    ctx.seq.resize(length);                                     // discard the tail so we can swap it into the memory
    ctx.generators_memory[(int16_t)ctx.periods.size()].swap(ctx.seq);
    ctx.seq.swap(ctx.seq_new);                                  // retrieve the new sequence from test_2
    ctx.seq.push_back((int16_t)ctx.c_cand);                     // add the candidates because we passed test_1 and test_2
    ctx.periods.push_back((int16_t)ctx.p_cand);
    int seq_size = (int)ctx.seq.size();
    ctx.seq_map[ctx.c_cand + length].push_back(seq_size - 1);
    int period = 0;
    while (true) {                                              // build the tail for this generator
        int curl = krul(ctx.seq, period, seq_size, 2);
        if (curl == 1)
            break;
        ctx.seq_map[curl + length].push_back(seq_size);
        ctx.seq.push_back((int16_t)curl);
        ++seq_size;
        ctx.periods.push_back((int16_t)period);
    }

    int tail = (int)ctx.periods.size();                         // find tail size (seq.size() - length)
    ctx.c_cand = 2;                                             // prepare candidates for next backtracking_step
    ctx.p_cand = 1 + tail / 2;
    ctx.change_indices.insert((int16_t)tail);                   // to improve generator, we would have to take note of the index where the 1 occured
    int len = real_generator_length(ctx);                       // retrieve actual generator length
    if (ctx.max_tails[len] < (int16_t)tail) {                   // and update maximum values for this thread
        ctx.max_tails[len] = (int16_t)tail;
        memcpy(&ctx.best_generators[len][0], &ctx.seq[length - len], sizeof(int16_t) * len);
    }
}

// check whether the current sequence allows for the candidates to be added
INLINING
bool test_1(context& ctx) {
    int l = (int)ctx.seq.size() - 1;                            // last element of pattern to check
    int lcp = l - ctx.p_cand;                                   // last element of potential pattern
    int limit = (ctx.c_cand - 1) * ctx.p_cand;                  // limit to which to check for repetition
    int16_t* p1 = &ctx.seq[l];
    int16_t* p2 = &ctx.seq[lcp];
    for (int i = 0; i < limit; ++i, --p1, --p2) {
        if (*p1 != *p2 and (*p1 | *p2) > 0)                     // check whether the repetition may be possible
            return false;
    }
    ctx.seq_new = ctx.seq;                                      // create dummy sequence
    ctx.pairs.clear();
    int pairs_size = 0;
    p1 = &ctx.seq_new[l];
    p2 = &ctx.seq_new[lcp];

    for (int i = 0; i < limit; ++i, --p1, --p2) {
        int16_t a = *p1;
        int16_t b = *p2;
        if (a != b) {                                           // because we are changing values below, we may encounter a possible break, again
            if ((a | b) > 0)
                return false;
            if (a > b)
                std::swap(a, b);                                // a is now always < b and < 0
            ctx.pairs.push_back(b);                             // add (b, a) combo to the pairs
            ctx.pairs.push_back(a);
            ++pairs_size;
            ++pairs_size;

            ctx.temp.clear();
            ctx.temp.push_back(a);                              // temporary vector that will hold all map values that need to be changed
            int temp_size = 1;
            for (int index = 0; index < temp_size; ++index) {   // because we don't (want to) change the map here (not sure we pass tests 1 and 2)
                int16_t* pi = &ctx.pairs[0];
                int16_t tmp = ctx.temp[index];
                for (int i = 0; i < pairs_size; i += 2, ++pi)   // if we change a to b, and later change b, we also need to change a in that case
                    if (*pi++ == tmp) {                         // so we need to check if we already crossed the value b
                        ctx.temp.push_back(*pi);
                        temp_size++;
                    }
            }
            for (int x : ctx.temp)                              // apply changes to the sequence
                for (int ind : ctx.seq_map[x + length])
                    ctx.seq_new[ind] = b;
        }
    }
    return true;
}

// check whether the proposed candidates invalidate the generator (regarding curl or period)
INLINING
bool test_2(context& ctx) {
    int l = (int)ctx.seq_new.size();
    int period = 0;

    for (int i = 0; i < l - length; ++i) {                                          // check within tail for a valid curl or period
        int curl = krul(ctx.seq_new, period, length + i, ctx.seq_new[length + i]);  // calculate curl and period up to this part of the sequence
        if (curl != ctx.seq_new[length + i] or period != ctx.periods[i])            // if the curl or its period are not related, the change is invalid
            return false;
    }
    int curl = krul(ctx.seq_new, period, l, ctx.c_cand);
    return (curl == ctx.c_cand and period == ctx.p_cand);                           // and also check the new candidate
}

// set a new step in the backtracking algorithm
INLINING
void backtracking_step(context& ctx) {
    if (test_1(ctx) && test_2(ctx)) // depending on whether the sequence will improve the tail length...
        append(ctx);                // we make the tail longer if it passed the checks
    else
        up(ctx);                    // or we upgrade the generator if they failed
}

int main(int argc, char *argv[])
{
    MPI_Init (&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // get rank of this process
	
    if (rank == 0) {        // master rank
        
        double t1_master = MPI_Wtime();
        FILE_OPEN;
        OUTPUT << "Length: " << length << std::endl;
        std::cout << "Hello from the master" << std::endl;
            
        int np;
        MPI_Comm_size(MPI_COMM_WORLD, &np); // get total number of processes
        
        int c1 = 2, p1 = 1, c2 = 2, p2 = 1, c3 = 2, p3 = 1;         // we support a depth of three (for additional fragmentation and thus scalability)
        int id;
        
        std::cout << "Distributing sequences" << std::endl << "2";
        // handle processes
        while(true) {
            int values[8];
            if (c1 > length) {                                              // end of program
                std::cout << std::endl << "Distributed all sequences!" << std::endl;
                break;
            }
            if ((c1 - 1) * p1 <= g_limit) {                                 // enter depth two
                if (c2 <= c1) {                                             // enter depth three
                    int new_values[8] = { c1, p1, 1, c2, p2, 1, c3, p3 };
                    memcpy(&values, new_values, sizeof(new_values));
                    if (++p3 > (length + 2) / c3) {
                        p3 = 1;
                        if (++c3 > length + 2) {
                            c3 = 2;
                            if (++p2 > (length + 1) / c2) {
                                p2 = 1;
                                c2++;
                            }
                        }
                    }
                }
                else {
                    int new_values[8] = { c1, p1, 1, c2, p2, 0, 0, 0 };
                    memcpy(&values, new_values, sizeof(new_values));
                    if (++p2 > (length + 1) / c2) {
                        p2 = 1;
                        if (++c2 > length + 1) {
                            c2 = 2;
                            if (++p1 > length / c1) {
                                p1 = 1;
                                c1++;
                                std::cout << std::endl << c1;
                            }
                            std::cout << " " << p1;
                        }
                    }
                }
            }
            else {
                int new_values[8] = { c1, p1, 0, 0, 0, 0, 0, 0 };
                memcpy(&values, new_values, sizeof(new_values));
                if (++p1 > length / c1) {
                    p1 = 1;
                    c1++;
                    std::cout << std::endl << c1;
                }
            }
            MPI_Recv(&id, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        // get notified of finished rank
            MPI_Send(&values[0], 8, MPI_INT, id, 0, MPI_COMM_WORLD);                                // send it new values
        }
        
        std::cout << "Gathering data";
        // clean up all processes one by one
        int values[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        int last_values[8];
        double elapsed;
        int max_tails[length + 1] = { 0 };
        int g_max_tails[length + 1] = { 0 };
        int16_t best_generators[length + 1][length];
        int16_t g_best_generators[length + 1][length];
        for (int i = 1; i < np; i++) {
            MPI_Recv(&id, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);         // get notified that this rank finished
            MPI_Send(&values[0], 8, MPI_INT, id, 0, MPI_COMM_WORLD);                                 // send it terminating values
            MPI_Recv(&elapsed, 1, MPI_DOUBLE, id, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);             // receive its elapsed time for logging
            MPI_Recv(&last_values[0], 8, MPI_INT, id, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            OUTPUT << "Rank: " << id << ", duration: " << elapsed;                      	         // log data
            for (int x : last_values)
                OUTPUT << " " << x;
            OUTPUT << std::endl;
            MPI_Recv(&max_tails[0], length + 1, MPI_INT, id, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive its maximum tails
            MPI_Recv(&(best_generators[0][0]), (length + 1)*length, MPI_INT16_T, id, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j <= length; ++j)                                                        // update global tails
                if (max_tails[j] > g_max_tails[j]) {
                    g_max_tails[j] = max_tails[j];
                    memcpy(&g_best_generators[j][0], &best_generators[j][0], sizeof(int16_t) * length);
                }
            std::cout << ".";
        }
        std::cout << std::endl;
            
        // process the tails
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
                    if (x != 0) 
                        OUTPUT << x << ",";
                OUTPUT << "]" << std::endl;
            }
        }
        std::cout << "Finished!" << std::endl;
        FILE_CLOSE;
        double t2_master = MPI_Wtime();
        elapsed = t2_master - t1_master;
        std::cout << "Took " << elapsed << " seconds." << std::endl;
        std::cout << "Quitting the program..." << std::endl;
	}
	
	else {                  // worker ranks
		
	    context ctx;
	    for (int i = 0; i <= length; ++i)                  	        // initiate local record vectors
            ctx.max_tails[i] = 0;
		
	    double t1_worker = MPI_Wtime();
		
	    int values[8];
        int last_values[8];
	    while (true) {
            MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);      // advertise that this rank is ready for a new combination
            MPI_Recv(&values[0], 8, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // get new values
            if (values[0] == 0)                                     // terminate if necessary
                break;
            int c1 = values[0], p1 = values[1];
            ctx.c2p2 = values[2];
            int c2 = values[3], p2 = values[4];
            ctx.c3p3 = values[5];

            ctx.change_indices.clear();
            
            for (int i = 0; i < length; ++i)                        // initiate sequence
                ctx.seq[i] = (int16_t)(-length + i);

            if (ctx.c2p2) {                                         // if we use depth 2, we add the candidates and change the sequence accordingly
                for (int i = 2; i <= c1; i++)
                    memcpy(&ctx.seq[length - i * p1], &ctx.seq[length - p1], p1 * sizeof(int16_t));
                ctx.seq.push_back(c1);
                ctx.periods.push_back(p1);
                ctx.change_indices.insert(0);
            }
            
            for (auto& v : ctx.seq_map)                             // clear the map
                v.clear();
            
            for (int j = 0; j < length + ctx.c2p2; ++j)             // and reconstruct
                ctx.seq_map[ctx.seq[j] + length].push_back(j);

            if (ctx.c3p3) {                                         // if we use depth 3, try to add the candidates after performing test_1 and test_2 for validity
                ctx.c_cand = c2;
                ctx.p_cand = p2;
                if (test_1(ctx) && test_2(ctx)) {
                    for (int i = 0; i < ctx.pairs.size(); i += 2) {             // NOW we are sure we passed test_1 and test_2, so it is time to update the map
                        for (int x : ctx.seq_map[ctx.pairs[i + 1] + length])
                            ctx.seq_map[ctx.pairs[i] + length].push_back(x);    // so we move all changed values to their new location
                        ctx.seq_map[ctx.pairs[i + 1] + length].clear();         // and delete their previous entries
                    }
                    ctx.seq.swap(ctx.seq_new);                                  // retrieve the new sequence from test_2
                    ctx.seq.push_back((int16_t)c2);                             // add the candidates because we passed test_1 and test_2
                    ctx.periods.push_back((int16_t)p2);
                    ctx.seq_map[c2 + length].push_back(length + 1);
                    ctx.change_indices.insert(1);
                }
                else {
                    ctx.seq.pop_back();
                    ctx.periods.pop_back();
                    continue;
                }
            }

            ctx.c_cand = values[(ctx.c2p2 + ctx.c3p3) * 3];         // select relevant candidates
            ctx.p_cand = values[(ctx.c2p2 + ctx.c3p3) * 3 + 1];

            backtracking_step(ctx);                                 // perform backtracking for this combination (c_cand, p_cand)
            while (ctx.periods.size() > ctx.c2p2 + ctx.c3p3)
                backtracking_step(ctx);

            if (ctx.c2p2) {                                         // clear the depth 2 and 3 candidates, if necessary
                ctx.seq.pop_back();
                ctx.periods.pop_back();
                if (ctx.c3p3) {
                    ctx.seq.pop_back();
                    ctx.periods.pop_back();
                }
            }
            memcpy(&last_values, values, sizeof(values));
        }
        double t2_worker = MPI_Wtime();
        double elapsed = t2_worker - t1_worker;
        MPI_Send(&elapsed, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);                // send elapsed time for logging
        MPI_Send(&last_values, 8, MPI_INT, 0, 2, MPI_COMM_WORLD);		        // send last values for debugging
        MPI_Send(&ctx.max_tails[0], length + 1, MPI_INT, 0, 3, MPI_COMM_WORLD); // send maximum tails in this rank
        MPI_Send(&(ctx.best_generators[0][0]), (length + 1)*length, MPI_INT16_T, 0, 4, MPI_COMM_WORLD);
    }   
    MPI_Finalize();
}