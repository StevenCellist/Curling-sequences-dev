// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 13-07-2021; completed time to length 48: 50 milliseconds*
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, GCC GNU compiler

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstring>
#include <algorithm>
#include <mpi.h>

#define INLINING inline __attribute__((always_inline))
typedef std::vector<int16_t> v16_t;

const int length = 120;     // Tweakable parameter: set this to the desired generator length (n)
const int limit = 30;       // Tweakable parameter: increase this value if ranks do not finish simultaneously (necessary for large # of ranks, preferable)
const int max_depth = 5;    // Tweakable parameter: controls fragmentation of variables; recommended value: largest power of 2 within max_np size
const int interval = length;// Tweakable parameter: set this to the desired interval between log ticks
const int max_np = 8;       // Tweakable parameter: set this to equal or more than the total number of processes in execution

struct context {                                            // all necessary variables for a rank;
    int c_cand = 0, p_cand = 0, depth = 1, max_tails[length + 1] = { 0 };
    v16_t seq = v16_t(length), seq_new, periods, pairs, temp;
    int16_t best_generators[length + 1][length] = { 0 };
    int generators_counter[length + 1] = { 0 };
    std::map<int16_t, v16_t> generators_memory;
    v16_t change_indices = { 0 };
    std::array<std::vector<int>, length> seq_map;
};

std::unordered_map<int, int> expected_tails = {             // the record tail lengths we found so far
    {2, 2},     {4, 4},     {6, 8},     {8, 58},    {9, 59},     {10, 60},    {11, 112},   {14, 118},   {19, 119},   {22, 120},
    {48, 131},  {68, 132},  {73, 133},  {77, 173},  {85, 179},   {115, 215},  {116, 228},  {118, 229},  {128, 332},  {132, 340},
    {133, 342}, {149, 343}, {154, 356}, {176, 406}, {197, 1668}, {199, 1669}, {200, 1670}, {208, 1708}, {217, 1836}, {290, 3382},
    {385, 3557}
};

// erase element from vector
void erase(std::vector<int>& v, int x) {
    int i = 0;
    while (v[i] != x)       // find requested element
        ++i;
    v[i] = v.back();        // move to last value
    v.pop_back();           // remove (now duplicated) last value
}

// specialized compare function using pointers
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

// find the shortest generator that could have generated current tail
int real_generator_length(context& ctx) {
    int i = 0;
    while ((ctx.seq[i] == (-length + i)) && (++i != length)) {}
    return (length - i);
}

INLINING // modify the candidate(s) and maybe the sequence to try and maximize the tail length
void up(context& ctx) {
    ++ctx.p_cand;                                                           // try period one larger now
    while (true) {
        if (ctx.periods.size() < ctx.depth)                                 // stop if tail is empty (relative to depth)
            break;
        if ((ctx.c_cand * ctx.p_cand) <= ctx.seq.size())                    // if this pair is within sequence size, we can break and try that instead
            break;
        if (ctx.c_cand > ctx.seq.back()) {                                  // if curl > last_curl, then raising curl with 1 would give curl > last_curl + 1.  
            int16_t k = (int16_t)ctx.periods.size() - 1;                    // due to a mathematical argument, this will give no solutions; therefore, the sequence will have to change to become longer.
            if (ctx.change_indices.back() == k + 1)                         // delete the current entry from change_indices if it was present
                ctx.change_indices.pop_back();    
            if (ctx.change_indices.back() == k) {                           // we will have to change the for-last curl in order to change the tail, so let's check whether it has been done yet                   
                ctx.c_cand = ctx.seq.back();                                // found it? let's try this same curl but now with one higher period
                ctx.p_cand = ctx.periods.back() + 1;
                v16_t& temp = ctx.generators_memory[k];                     // retrieve generator from memory
                for (int i = 0; i < length; i++)
                    if (ctx.seq[i] != temp[i]) {                            // check where the differences are, apply them, and change map accordingly
                        if (ctx.seq[i] < 0)
                            erase(ctx.seq_map[ctx.seq[i] + length], i);
                        ctx.seq_map[temp[i] + length].push_back(i);
                        ctx.seq[i] = temp[i];
                    }
                ctx.generators_memory.erase(k);                             // and delete memory entry
            }
            else {                                                          // didn't find it? then using a mathematical argument, we can skip a number of cases. 
                ctx.change_indices.push_back(k);                            // insert that index and let's try increased curl with matching period
                ctx.c_cand = ctx.seq.back() + 1;
                ctx.p_cand = 1 + k / ctx.c_cand;
            }
            ctx.seq.pop_back();
            ctx.periods.pop_back();
        }
        else {
            ++ctx.c_cand;                                                   // if not, we increase the curl to try and calculate its period
            ctx.p_cand = 1 + (int)ctx.periods.size() / ctx.c_cand;
        }
    }
}

INLINING // construct the tail for this generator
void append(context& ctx) {
    for (int i = 0; i < ctx.pairs.size(); i += 2) {             // now we are sure we passed test_1 and test_2, so it is time to update the map
        if (ctx.pairs[i] < 0)
            for (int x : ctx.seq_map[ctx.pairs[i + 1] + length])
                ctx.seq_map[ctx.pairs[i] + length].push_back(x);// so we move all changed values to their new location
        ctx.seq_map[ctx.pairs[i + 1] + length].clear();         // and delete their previous entries
    }
    ctx.seq.resize(length);                                     // discard the tail so we can swap it into the memory
    ctx.generators_memory[(int16_t)ctx.periods.size()].swap(ctx.seq);
    ctx.seq.swap(ctx.seq_new);                                  // retrieve the new sequence from test_2
    ctx.seq.push_back((int16_t)ctx.c_cand);                     // add the candidates because we passed test_1 and test_2
    ctx.periods.push_back((int16_t)ctx.p_cand);
    int seq_size = (int)ctx.seq.size();
    int period = 0;
    while (true) {                                              // build the tail for this generator
        int curl = krul(ctx.seq, period, seq_size, 2);
        if (curl == 1)
            break;
        ctx.seq.push_back((int16_t)curl);
        ++seq_size;
        ctx.periods.push_back((int16_t)period);
    }

    int tail = (int)ctx.periods.size();                         // find tail size (seq.size() - length)
    ctx.c_cand = 2;                                             // prepare candidates for next backtracking_step
    ctx.p_cand = 1 + tail / 2;
    ctx.change_indices.push_back((int16_t)tail);                // to improve generator, we would have to take note of the index where the 1 occured
    int len = real_generator_length(ctx);                       // retrieve actual generator length
    if (ctx.max_tails[len] < (int16_t)tail) {                   // and update maximum values for this thread
        ctx.generators_counter[len] = 1;                        // (re)set counter for this record
        ctx.max_tails[len] = (int16_t)tail;
        memcpy(&ctx.best_generators[len][0], &ctx.seq[length - len], sizeof(int16_t) * len);
    }
    else if (ctx.max_tails[len] == (int16_t)tail)
        ctx.generators_counter[len]++;                          // increase counter for this record
}

INLINING // check whether the current sequence allows for the candidates to be added
bool test_1(context& ctx) {
    int l = (int)ctx.seq.size() - 1;                            // last element of pattern to check
    int lcp = l - ctx.p_cand;                                   // last element of potential pattern
    int limit = (ctx.c_cand - 1) * ctx.p_cand;                  // limit to which to check for repetition
    int16_t* p1 = &ctx.seq[l];
    int16_t* p2 = &ctx.seq[lcp];
    for (int i = 0; i < limit; ++i, --p1, --p2)
        if (*p1 != *p2 and (*p1 | *p2) > 0)                     // check whether the repetition may be possible
            return false;
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
            pairs_size += 2;

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

INLINING // check whether the proposed candidates invalidate the generator (regarding curl or period)
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

INLINING // set a new step in the backtracking algorithm
void backtracking_step(context& ctx) {
    if (test_1(ctx) && test_2(ctx))     // depending on whether the sequence will improve the tail length...
        append(ctx);                    // we make the tail longer if it passed the checks
    else
        up(ctx);                        // or we upgrade the generator if they failed
}

// main rank
void master(const int rank, const int np) {
    double t1 = MPI_Wtime();
    int id, values[1 + 2 * max_depth];                                              // initiate with default values
    values[0] = max_depth;
    for (int i = 1; i <= max_depth; i++) {
        values[2 * i - 1] = 2;
        values[2 * i] = 1;
    }
    int depth = max_depth;
    std::cout << "Master: distributing sequences...\n2";
    while (true) {
        if (values[depth * 2 - 1] * values[depth * 2] < length + depth) {           // check if candidates are within sequence size
            MPI_Recv(&id, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    // get notified of finished worker
            MPI_Send(&id, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);                                    // send that ID to the logger
            MPI_Send(&values[0], 1 + 2 * max_depth, MPI_INT, 1, 2, MPI_COMM_WORLD);             // as well as the values
            MPI_Send(&values[0], 1 + 2 * max_depth, MPI_INT, id, 2, MPI_COMM_WORLD);            // send that worker new values
            values[depth * 2]++;                                                    // and increase period candidate by one
        }
        else {
            values[depth * 2] = 1;                                              // if this combination of curl and period was too large, reset period
            values[depth * 2 - 1]++;                                            // and increase curl by one
            if (depth > 1) {
                if (values[depth * 2 - 1] > values[depth * 2 - 3] + 1) {        // curl larger than previous curl + 1? skip to next (mathematical)   
                    values[depth * 2 - 1] = 2;                                  // reset the curl of the current depth
                    values[depth * 2 - 2]++;                                    // and increase the period of previous depth with one
                    if (depth == 2) 
                        std::cout << " " << values[2];
                }
            }
            else {
                if (values[1] > length)                                             // curl too large for sequence size?
                    break;                                                          // then we just gave out the last set of values
                if (values[1] <= limit) 
                    std::cout << "\b\b" << "  " << std::endl << values[1];
            }
            int sum = values[1] * values[2];                                        // depth is at least one to get us started on the sum
            for (depth = 1; depth < max_depth; depth++) {                           // check the (weighed) depth that matches current values
                if (sum > limit)
                    break;                                                          // if we go over the limit, we do not do this combination
                int next = (depth + 1) * values[depth * 2 + 1] * values[depth * 2 + 2];
                sum += next;                                                        // otherwise, we accept this depth
            }
            values[0] = depth;                                                      // store the current depth
        }
    }
    std::cout << "\b\b" << "  " << "\nMaster: terminating loggers and workers. " << MPI_Wtime() - t1 << " seconds." << std::endl;
    values[0] = 0;
    for (int i = 3; i < np; i++) {
        MPI_Recv(&id, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    // get notified that a worker has finished
        MPI_Send(&values[0], 1 + 2 * max_depth, MPI_INT, id, 2, MPI_COMM_WORLD);            // send it terminating values
    }
    MPI_Send(&rank, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);                              // terminate value logger
    MPI_Recv(&id, 1, MPI_INT, 1, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            // wait for value logger termination
    MPI_Recv(&id, 1, MPI_INT, 2, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            // wait for results logger termination
    std::cout << MPI_Wtime() - t1 << " seconds.\n";
    std::cout << "Master: finished.\n";                                             // terminate master
}

// log the rank values to file upon interval tick
void value_logger(const int rank, const int np) {
    std::ofstream log_file;
    int cycles = np * interval;
    int id, log[max_np][1 + 2 * max_depth], values[1 + 2 * max_depth];
    while (true) {
        MPI_Recv(&id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                         // receive ID from master
        if (id == 0)    // termination signal from master
            break;
        MPI_Recv(&values[0], 1 + 2 * max_depth, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // receive values from master
        memcpy(&log[id][0], &values[0], sizeof(values));                                            // store values into log
        cycles--;
        if (cycles == 0) {                                  // if we went through 'interval' cycles, write the log to file
            log_file.open("Log.txt");
            log_file << "Length: " << length << std::endl;
            for (int i = 3; i < np; i++) {
                log_file << i << ", depth: " << log[i][0] << " = ";
                for (int j = 1; j <= log[i][0]; j++)
                    log_file << "(" << log[i][j * 2 - 1] << ", " << log[i][j * 2] << ") ";
                log_file << std::endl;
            }
            log_file.close();
            cycles = np * interval;
        }
    }
    MPI_Send(&rank, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);     // report termination to master
    std::cout << "Value logger: finished.\n";               // terminate logger
}

INLINING // write relevant results to file
void write_results(int (&g_max_tails)[length + 1], int16_t (&g_best_generators)[length + 1][length], int (&g_generators_counter)[length + 1]) {
    std::ofstream results_file;
    results_file.open("Results.txt");
    results_file << "Length: " << length << std::endl;
    int record = 0;
    for (int i = 0; i <= length; ++i)
        if (g_max_tails[i] > record) {
            record = g_max_tails[i];
            if (expected_tails.find(i) == expected_tails.end())
                results_file << "NEW:" << std::endl;
            else if (expected_tails[i] != record)
                results_file << "WRONG:" << std::endl;
            results_file << i << ": " << record << ", [";
            for (int x : g_best_generators[i])
                if (x != 0) 
                    results_file << x << ",";
            results_file << "]" << std::endl;
            if (g_generators_counter[i] > 1)
                results_file << "Found non-unique generator: n = " << i << " with frequency " << g_generators_counter[i] << "!\n";
        }
    results_file.close();
}

// log the results to file upon interval tick and termination
void results_logger(const int rank, const int np) {
    double t1 = MPI_Wtime();
    bool running = true;
    int id = 0, running_processes = np - 3;
    int last_values[1 * 2 * max_depth], max_tails[length + 1] = { 0 }, g_max_tails[length + 1] = { 0 }, generators_counter[length + 1] = { 0 }, g_generators_counter[length + 1] = { 0 };
    int16_t best_generators[length + 1][length], g_best_generators[length + 1][length];
    std::ofstream ranks_file;
    ranks_file.open("Ranks.txt");
    ranks_file << "Length: " << length << std::endl;
    while (running_processes) {
        MPI_Recv(&id, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                                // get notified that a worker is terminating
        MPI_Recv(&running, 1, MPI_C_BOOL, id, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_tails[0], length + 1, MPI_INT, id, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                         // receive its maximum tails
        MPI_Recv(&(best_generators[0][0]), (length + 1)*length, MPI_INT16_T, id, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // as well as corresponding generators
        MPI_Recv(&generators_counter[0], length + 1, MPI_INT, id, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                // and counters for the generators
        for (int k = 0; k <= length; k++) {                                                                             // update current best tails
            if (max_tails[k] > g_max_tails[k]) {
                if (!running) g_generators_counter[k] = generators_counter[k];
                g_max_tails[k] = max_tails[k];
                memcpy(&g_best_generators[k][0], &best_generators[k][0], sizeof(int16_t) * length);
            }
            else if (max_tails[k] == g_max_tails[k])
                if (!running) g_generators_counter[k] += generators_counter[k];                                         // increase counter if we found non-unique generators
        }
        MPI_Recv(&last_values[0], 1 + 2 * max_depth, MPI_INT, id, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                // receive its last values for logging
        if (!running) {
            running_processes--;
            ranks_file << "Rank: " << id << ", duration: " << (MPI_Wtime() - t1) << ", " << last_values[0] << ": ";     // log worker number, duration and last values
            for (int j = 1; j <= last_values[0]; j++) 
                ranks_file << "(" << last_values[j * 2 - 1] << ", " << last_values[j * 2] << "), ";
            ranks_file << std::endl;
            std::cout << "*";
        }
        else
            write_results(g_max_tails, g_best_generators, g_generators_counter);          // write final results to file
    }
    ranks_file.close();
    write_results(g_max_tails, g_best_generators, g_generators_counter);          // write final results to file
    MPI_Send(&rank, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);     // report termination to master
    std::cout << "\nResults logger: finished.\n";           // terminate logger
}

void worker(const int rank, const int np) {
    context ctx;
    int cycles = interval;
    int values[1 + 2 * max_depth], last_values[1 + 2 * max_depth];
    while (true) {
        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);                  // advertise that this rank is ready for a new combination
        MPI_Recv(&values[0], 1 + 2 * max_depth, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // get new values
        bool running = values[0];
        if (!running or !cycles) {
            MPI_Send(&rank, 1, MPI_INT, 2, 3, MPI_COMM_WORLD);                                              // send rank number to results logger
            MPI_Send(&running, 1, MPI_C_BOOL, 2, 4, MPI_COMM_WORLD);                                        // send whether this worker is still going
            MPI_Send(&ctx.max_tails[0], length + 1, MPI_INT, 2, 5, MPI_COMM_WORLD);                         // send maximum tails of this worker
            MPI_Send(&(ctx.best_generators[0][0]), (length + 1)*length, MPI_INT16_T, 2, 6, MPI_COMM_WORLD); // send best generators of this worker
            MPI_Send(&ctx.generators_counter[0], length + 1, MPI_INT, 2, 7, MPI_COMM_WORLD);                // send counters for found records
            MPI_Send(&last_values, 1 + 2 * max_depth, MPI_INT, 2, 8, MPI_COMM_WORLD);                       // send last values of this worker
            if (!running)
                break;
            cycles = interval;
        }
        cycles--;
        ctx.depth = values[0];                                              // store depth as variable

        if (last_values[0] > 1) {
            ctx.seq.resize(length);                                         // reset to default values
            ctx.periods.clear();
            ctx.change_indices.clear();
        }
        memcpy(&last_values, values, sizeof(values));                       // store current values for logging

        for (int i = 0; i < length; ++i)                                    // initiate sequence
            ctx.seq[i] = (int16_t)(-length + i);

        if (ctx.depth > 1) {                                                // if we use depth 2 or higher, we add the first candidates and change the sequence accordingly
            for (int i = 2; i <= values[1]; i++)
                memcpy(&ctx.seq[length - i * values[2]], &ctx.seq[length - values[2]], values[2] * sizeof(int16_t));
            ctx.seq.push_back((int16_t)values[1]);
            ctx.periods.push_back((int16_t)values[2]);                      // we don't change generators_memory or change_indices, since we won't need them for these low values of i
        }

        for (auto& v : ctx.seq_map)                                         // clear the map
            v.clear();
        for (int j = 0; j < length; ++j)                                    // and reconstruct
            ctx.seq_map[ctx.seq[j] + length].push_back(j);

        bool next = false;
        for (int i = 3; i <= ctx.depth; i++) {                              // if we use depth 3 or more, try to add the candidates after performing test_1 and test_2 for validity
            ctx.c_cand = values[i * 2 - 3];
            ctx.p_cand = values[i * 2 - 2];
            if (test_1(ctx) && test_2(ctx)) {
                for (int i = 0; i < ctx.pairs.size(); i += 2) {             // now we are sure we passed test_1 and test_2, so it is time to update the map
                    if (ctx.pairs[i] < 0) {
                        for (int x : ctx.seq_map[ctx.pairs[i + 1] + length])
                            ctx.seq_map[ctx.pairs[i] + length].push_back(x);    // so we move all changed values to their new location
                    }
                    ctx.seq_map[ctx.pairs[i + 1] + length].clear();         // and delete their previous entries
                }
                ctx.seq.swap(ctx.seq_new);                                  // retrieve the new sequence from test_2
                ctx.seq.push_back((int16_t)ctx.c_cand);                     // add the candidates because we passed test_1 and test_2
                ctx.periods.push_back((int16_t)ctx.p_cand);                 // we don't change generators_memory or change_indices, since we won't need them for these low values of i
            }
            else {
                next = true;                                                // if some combination failed, skip and request new variables
                break;
            }
        }
        if (next)
            continue;
        ctx.c_cand = values[ctx.depth * 2 - 1];                             // select relevant candidates
        ctx.p_cand = values[ctx.depth * 2];

        do backtracking_step(ctx); 
        while (ctx.periods.size() >= ctx.depth);                            // perform backtracking for this combination (c_cand, p_cand)
    }
}

int main(int argc, char *argv[])
{
    int rank, np;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // get rank number of this process
    MPI_Comm_size(MPI_COMM_WORLD, &np);     // get total number of processes
    
    if (rank == 0)      // master rank
        master(rank, np);
    else if (rank == 1) // values logger
        value_logger(rank, np);
    else if (rank == 2) // results logger
        results_logger(rank, np);
    else                // worker ranks
        worker(rank, np);  
    MPI_Finalize();
}