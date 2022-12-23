// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 07-06-2021; completed time to length 48: 68 milliseconds*
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, Microsoft VS Studio 2019 Compiler

#define CPPHTTPLIB_OPENSSL_SUPPORT
#include "httplib.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <array>
#include <thread>
#include <mutex>
#include <queue>

#ifdef _MSC_VER
#define FILE 
#define FILE_OPEN   
#define OUTPUT stdout
#define FILE_CLOSE  
#define INLINING __declspec(noinline)
#define UPDATE_TIME ctime_s(date, 26, &now)
#else // gcc (Linux)
#include <cstring>
#include <algorithm>
#define FILE // std::ofstream file
#define FILE_OPEN // file.open("Sequences.txt")
#define OUTPUT stdout // file
#define FILE_CLOSE // file.close()
#define INLINING __attribute__((always_inline))
#define UPDATE_TIME ctime_r(&now, date)
#endif

using namespace std::chrono;
typedef std::vector<int16_t> v16_t;

/*
    @brief: all necessary variables required for a worker to construct curling sequences using backtracking
*/
struct context {
    int length;		// user-selected length. The maximal length of the generator.
    int c_cand;		// next curling number candidate to append to sequence
    int p_cand;		// next period candidate for the curling number candidate
    int depth;		// depth used for thread-splitting
    v16_t seq;		// curling sequence
    v16_t seq_new;	// curling sequence that may or may not be used
    v16_t periods;	// periods for all curling numbers in the tail
    v16_t pairs;	// pairs of numbers in the sequence that may or may not need to be updated in seq_map (may hold duplicates)
    v16_t temp;		// pairs of numbers in the sequence that may or may not need to be updated in seq_map (no duplicates)
    std::vector<v16_t> seq_map;		            // reverse map of seq (value + ctx.length -> [indices] instead of index -> value)
                                                // we view the elements of seq_map as sets, but in this code we use vectors for efficiency  
    std::unordered_set<int16_t> change_indices;	// all indices where we branch off in backtracking as we will return here later
    std::map<int16_t, v16_t> grts_mem;	        // all matching generators at points where we branch off in backtracking as we will return there later
    std::vector<int> best_tails;				// all longest tail lengths (per worker)
    std::vector<v16_t> best_grts;			    // all generators for the longest tail lengths (per worker)
};

/*
    @brief: all global variables in order to synchronize between workers
*/
std::queue<v16_t> cmbs;					        // queue with the combinations to process by workers

/*
    @brief: compare two chunks of memory
    in order to compare two chunks of memory as fast as possible, create pointers to 64 bits at once - 4x int16_t
    after that, compare any residue (16, 32 or 48 bits) using 'normal' int16_t pointers

    input:
    const int16_t* p1 - pointer to first 16 bits of memory block 1
    const int16_t* p2 - pointer to first 16 bits of memory block 2
    int count         - number of values to compare (int16_t)

    returns:
    bool diff - true in case of any difference / false in case of 100% match
*/
bool diff(const int16_t* p1, const int16_t* p2, int count) {
    int count64 = count / 4;                    // number of four x 16-bit ints in uint64_t
    if (count64) {
        uint64_t* p1_64 = (uint64_t*)p1;        // cast to 64-bit pointers
        uint64_t* p2_64 = (uint64_t*)p2;
        while (count64--)                       // loop through all 64-bit blocks
            if (*p1_64++ != *p2_64++)
                return true;                    // found a difference
        p1 = (int16_t*)p1_64;                   // cast back to 16-bit pointer for residue
        p2 = (int16_t*)p2_64;
    }

    int count16 = count % 4;                    // residue
    if (count16) {
        while (count16--)                       // loop through last 16-bit blocks
            if (*p1++ != *p2++)
                return true;                    // found a difference
    }
    return false;                               // found no single difference
}

/*
    @brief: find the highest frequency of a sequence following the curling principle
    the function 'krul' (Dutch for 'curl') is a slight generalization of the curling number function.
    for a string 's' and parameters 'len' and 'minimum', it returns max(cn(suffix),minimum-1).
    here 'suffix' is the suffix of 'seq' with length 'len', and cn(suffix) is the curling number of 'suffix'.
    if 'len' equals the length of 'seq' and if 'minimum' equals 2, then this function calculates the standard curling number of 'seq'.
    if the resulting 'curl' is at least 'minimum', 'period' is set to the smallest period corresponding to a 'curl' amount of repetitions at the end of 'seq'.

    input:
    const v16_t& seq - sequence to check for repetition
    int& period - integer to save the period
    int len - the length of the sequence
    int minimum - the minimum curl we are sure to find (only used when backtracking)

    returns:
    int curl - the [maximum of (minimum - 1) and] curling number of the prefix of seq with length len
    (indirectly) int& period - the minimal period corresponding to this curling number
*/
int krul(const v16_t& seq, int& period, int len, int minimum) {
    int curl = minimum - 1;                                 // base value for curl (only used when backtracking - otherwise default 1)
    int limit = len / minimum;                              // maximum period size
    const int16_t* p1 = &seq[len - 1];                      // pointer to the start of the last pattern
    for (int i = 1; i <= limit; ++i, --p1) {                // check for repetition up to the limit, starting with period = 1
        const int16_t* p2 = p1 - i;                         // pointer to start of the previous pattern
        for (int freq = 2; p2 >= &seq[0]; ++freq, p2 -= i) {
            if (diff(p1, p2, i))
                break;                                      // if it doesn't match, continue to next period
            if (curl < freq) {                              // found better curl?
                curl = freq;                                // update curl
                limit = len / (curl + 1);                   // update limit
                period = i;                                 // update period
            }
        }
    }
    return curl;                                            // return largest found curl
}


/*
    @brief: this function checks whether the current sequence allows for the candidates to be added.
    so we want to check whether c_cand and p_cand can be added.
    if this happens, then the last 'c_cand' blocks of length 'p_cand' of ctx.seq have to be equal.
    for this, we have to change the sequence: certain negative elements will have to be set to less-negative elements or positive elements
    to do this in a systematic way, we essentially create equivalence classes:
    two numbers are in the same equivalence class if they 'should be equal' according to the new repetition that we want to create.
    however, when two different positive numbers are placed in the same equivalence class, we obtain a contradiction, so we return false.
    if everything goes well, every equivalence class contains at most 1 positive number.
    every number in ctx.seq is then replaced by the highest element of its equivalence class.
    the resulting seqeunce is stored in ctx.seq_new (if any later test fails, we can continue from ctx.seq and discard ctx.seq_new)

    the implementation may look a bit weird, but it is fast

    input:
    context& ctx - all relevant variables
*/
INLINING
bool test_cands(context& ctx) {
    int l = (int)ctx.seq.size() - 1;                            // l and lcp are the last elements of the two blocks we need to check
    int lcp = l - ctx.p_cand;
    int limit = (ctx.c_cand - 1) * ctx.p_cand;                  // number of checks we have to make
    int16_t* p1 = &ctx.seq[l];                                  // fast pointer to indices
    int16_t* p2 = &ctx.seq[lcp];
    for (int i = 0; i < limit; ++i, --p1, --p2)                 // this initial, incomplete check already filters out a lot of the cases where the candidates are impossible
        if (*p1 != *p2 && (*p1 | *p2) > 0)                      // if two distinct positive numbers should be equal, we reach a contradiction
            return false;
    ctx.seq_new = ctx.seq;                                      // having completed the initial test, we create a dummy sequence to which we will apply changes.
    ctx.pairs.clear();                                          // below, we will need a map for ctx.seq_new.
                                                                // however, creating this new map costs time.
                                                                // we use a faster method instead, which is ctx.pairs.
                                                                // ctx.pairs is an auxiliary vector which will contain all pairs (b,a) where a < b and a has to be changed to b.
    ctx.pairs.reserve(limit);
    int pairs_size = 0;
    p1 = &ctx.seq_new[l];                                       // reset pointer to last elements
    p2 = &ctx.seq_new[lcp];

    for (int i = 0; i < limit; ++i, --p1, --p2) {
        int16_t a = *p1;                                        // the values of a and b should become equal
        int16_t b = *p2;
        if (a != b) {                                           // because we are changing values below, we may encounter a possible break, again.
            if ((a | b) > 0)                                    // if a and b are positive and distinct, we reach a contradiction
                return false;
            if (a > b)
                std::swap(a, b);                                // a is now always < b and < 0
            ctx.pairs.push_back(b);                             // add (b, a) combo to ctx.pairs
            ctx.pairs.push_back(a);
            pairs_size += 2;

            int16_t* p_begin = &ctx.pairs[0];
            int16_t* p_end = p_begin + pairs_size;

            ctx.temp.clear();                                   // temporary vector that will hold all numbers that are already changed to a
            ctx.temp.push_back(a);
            int temp_size = 1;
            for (int index = 0; index < temp_size; ++index) {   // because we don't (want to) change the map here (not sure we pass test_cands and test_seq_new).
                                                                // the goal of this for-loop is to find all values which have already been changed to a in seq_new.
                                                                // these values will be stored in ctx.temp.
                                                                // for each number in temp, we will loop through ctx.pairs to find all values which have been changed to that number;
                                                                // those values are then also added to temp.
                int16_t* pi = p_begin;
                int16_t tmp = ctx.temp[index];
                for (; pi < p_end; ++pi)                        // we loop through ctx.pairs
                                                                // note that 'pi' is raised twice in each iteration
                                                                // that is because we only want to see the first element of each pair
                    if (*(pi++) == tmp) {
                        ctx.temp.push_back(*pi);                // note that the element that we add to ctx.temp is the second element of the pair ('pi' has been raised in the line above)
                        temp_size++;
                    }
            }
            for (int x : ctx.temp)                              // change every 'a' in ctx.seq_new to a 'b'
                                                                // note that ctx.seq_map does correspond to ctx.seq, not to ctx.seq_new
                                                                // therefore we have to loop through ctx.temp
                for (int ind : ctx.seq_map[x + ctx.length])
                    ctx.seq_new[ind] = b;
        }
    }
    return true;                                                // if we did not reach a contradiction, then we continue to test_seq_new.
}

/*
    @brief: this function checks whether the proposed changed sequence ctx.seq_new does still work (regarding its curl and periods)
    of course, we only have to check this for elements of the tail.

    input:
    context& ctx - all relevant variables
*/
INLINING 
bool test_seq_new(context& ctx) {
    int l = (int)ctx.seq_new.size();
    int period = 0;
    for (int i = 0; i < l - ctx.length; ++i) {                                              // check within tail whether the curls and periods are valid
        int curl = krul(ctx.seq_new, period, ctx.length + i, ctx.seq_new[ctx.length + i]);  // calculate curl and minimal corresponding period up to this part of the sequence
        auto index = ctx.change_indices.find(i);                                            // check whether the sequence has been changed when this number was added
                                                                                            // this has to do with a complicated mathematical optimalization
                                                                                            // we explain this in our article
                                                                                            // TODO explain this in the article
        if (index == ctx.change_indices.end()) {                                            // didn't find it?
            if (curl != ctx.seq_new[ctx.length + i])                                        // if the curl is not equal to the values in ctx.seq_new, the change is invalid
                                                                                            // notice that the function krul actually calculates a maximum of the curling number and the expected number minus 1
                                                                                            // this is faster and still works in order to check whether the values are equal
                return false;
        }
        else {                                                                              // if the sequence has changed at this point,
                                                                                            // then our complicated optimalization tells us that we can make an extra requirement
            if (curl != ctx.seq_new[ctx.length + i] | period != ctx.periods[i])             // now both the curl and the minimal period have to be equal to the values in ctx.seq_new, in order to pass the test
                return false;
        }
    }
    int curl = krul(ctx.seq_new, period, l, ctx.c_cand);
    return (curl == ctx.c_cand && period == ctx.p_cand);                                    // and also check the new candidates
                                                                                            // since the next location will always end up in change_indices,
                                                                                            // we can use the complicated optimalization here also    
}

/*	main-thread function
    @brief: this function generates curl/period combinations for all threads to work on
    the combinations are put into a queue of size length*length to limit memory consumption
    every 100 milliseconds, it will fill the queue up to this size again

    a combination has the following structure:
    [depth, c_1, p_1, c_2, p_2, ..., c_n, p_n]
    where n = max_depth

    each combination is checked to make sure it produces even a single result
    if it is not a valid combination (i.e. len(tail) = 0), don't bother putting it in the queue
    (number of valid combinations is of order 1%)

    input:
    int len - the user-selected length
    int max_depth - the calculated maximum depth for this length
*/
void generate_combinations(int len, int max_depth, int n_lines) {

    httplib::Client cli("https://139027.icvi.nl");
    cli.enable_server_certificate_verification(false);
    std::string suffix = "/cs_insert.php/?length=";
    
    // send request to reset all results up to and including current length
    auto res = cli.Get(suffix += std::to_string(len)); 

    char date[26] = {};
    time_t now = time(NULL);
    UPDATE_TIME;

    int tries = 0;                              // retry counter
    while (res->body != "true") {               // retry until "true" is replied
        // std::cout << res->body << std::endl;    // if it failed, print error message 
        tries++;                                // increase retry counter
        if (tries > 10) {                       // if we tried more than 10 times, just give up
            fprintf(OUTPUT, "[%.15s] Timed out resetting results, continuing\n", &date[4]);
            break;
        }
        res = cli.Get(suffix += std::to_string(len));   // resend same request

        now = time(NULL);
        UPDATE_TIME;
        fprintf(OUTPUT, "[%.15s] Failed resetting results, retrying in %d seconds\n", &date[4], tries);
        std::this_thread::sleep_for(seconds(tries));    // incremental sleep
    }

    /*
        initialize all variables that are required for constructing the curling sequences
        resize all variables that have fixed or known minimum size to prevent unnecessary push_backs
    */
    context ctx;
    ctx.length = len;					// save length
    ctx.seq.resize(len);				// minimum sequence size is length
    ctx.seq_map.resize(2 * len + 2);    // adjust for +1 index on positive and negative side to simplify accessing elements
    ctx.best_tails.resize(len + 1);		// adjust for +1 index to simplify accessing elements
    ctx.best_grts.resize(len + 1);		// adjust for +1 index to simplify accessing elements
    for (int i = 0; i <= len; i++)
        ctx.best_grts[i].resize(len);	// generator is at most length 'len'

    /*
        the starting depth is the maximum depth
        and all curl/period combinations are at least 2/1
    */
    uint16_t depth = max_depth;
    v16_t cmb(1 + 2 * max_depth);
    cmb[0] = max_depth;
    for (int i = 1; i <= max_depth; i++) {
        cmb[2 * i - 1] = 2;
        cmb[2 * i] = 1;
    }

    // fill up the queue with all combinations until the first curl (depth 1) exceeds the sequence length (invalid)
    while (cmb[1] <= ctx.length) {
        if (cmbs.size() == n_lines) {
            // send HTTPS GET
            std::string tail = "";
            for (int j : cmbs.front()) {
                tail += std::to_string(j) + " ";
            }
            tail.pop_back();
            std::string suffix = "/cs_insert.php/?len=" + std::to_string(len) + "&num=" + std::to_string(n_lines) + "&combination=";
            auto res = cli.Get(suffix += tail);

            tries = 0;
            while (res->body != "true") {
                std::cout << res->body << std::endl;
                tries++;
                if (tries > 10) {
                    fprintf(OUTPUT, "[%.15s] Timed out uploading job, continuing\n", &date[4]);
                    std::cout << tail << std::endl;
                    break;
                }
                res = cli.Get(suffix += tail);

                now = time(NULL);
                UPDATE_TIME;
                fprintf(OUTPUT, "[%.15s] Failed uploading job, retrying in 1 second\n", &date[4]);
                std::this_thread::sleep_for(seconds(1));
            }

            std::queue<v16_t>().swap(cmbs);
        }
        else {
            /*
                reset all variables from previous calculations (much faster than destroying and reconstructing!)
                and then check if the current combination is valid
                in case it is valid, we store it in the queue
                qualitative explanation is given in the backtracking() function for the workers - here, they are just helper functions
            */
            ctx.depth = depth;                                          // set to current value
            ctx.seq.resize(ctx.length);
            for (int i = 0; i < ctx.length; ++i)
                ctx.seq[i] = (int16_t)(-ctx.length + i);
            for (auto& v : ctx.seq_map)
                v.clear();
            for (int j = 0; j < ctx.length; ++j)
                ctx.seq_map[j].push_back((int16_t)j);
            ctx.periods.clear();
            ctx.change_indices.clear();

            bool invalid = false;										// fast skip in case a combination is invalid
            // for all 'intermediate' depths, check if its combination yields a result 
            // and update the necessary variables since they influence the next depth-combination
            for (int i = 1; i < ctx.depth; i++) {
                ctx.c_cand = cmb[i * 2 - 1];							// load curl/period combination to struct
                ctx.p_cand = cmb[i * 2];
                if (test_cands(ctx) && test_seq_new(ctx)) {			    // check if we yield any result (else invalid combination)
                    for (int j = 0; j < ctx.pairs.size(); j += 2) {     // if it's valid, update necessary variables
                        for (int x : ctx.seq_map[ctx.pairs[j + 1] + ctx.length])
                            ctx.seq_map[ctx.pairs[j] + ctx.length].push_back(x);
                        ctx.seq_map[ctx.pairs[j + 1] + ctx.length].clear();
                    }
                    ctx.seq.swap(ctx.seq_new);
                    ctx.seq.push_back((int16_t)ctx.c_cand);
                    ctx.periods.push_back((int16_t)ctx.p_cand);
                    ctx.seq_map[ctx.c_cand + ctx.length].push_back((int16_t)(ctx.length + 1));
                    ctx.change_indices.insert((int16_t)(i - 1));
                }
                else {
                    invalid = true;                                     // if some combination yielded no result, we can skip to the next one
                    break;
                }
            }

            // for the final depth, we only need to check if it yields a result
            // as we won't store its matching variables: they are reconstructed by the workers
            ctx.c_cand = cmb[ctx.depth * 2 - 1];
            ctx.p_cand = cmb[ctx.depth * 2];
            if (!invalid && test_cands(ctx) && test_seq_new(ctx)) {

                // if this combination is valid, store its depth, lock the queue to avoid conflicts, and push it to the back of the queue
                cmb[0] = depth;
                cmbs.push(cmb);
            }

            /*
                continue to the next combination
                the default behaviour is to increase the period of the current depth by one
                but, if c_depth x p_depth > len(sequence), we go out of bounds and must reset the period and increase the curl by one
                but, if c_depth > c_{depth - 1}, we can profit from a mathematical shortcut: this will by definition yield zero results
                in that case, we go one up in depth and increase its period by one (starting the same reasoning again)
            */
            cmb[depth * 2]++;                                          	        // increase period by one
            bool recalc = false;												// keep track whether we need to recalculate the depth
            while (cmb[depth * 2 - 1] * cmb[depth * 2] >= ctx.length + depth) {	// check whether outside sequence size
                cmb[depth * 2] = 1;                                             // reset current period to default
                cmb[depth * 2 - 1]++;                                           // and increase current curl by one
                recalc = true;													// we will need to recalculate depth after this since we lowered a value
                if (depth > 1 && cmb[depth * 2 - 1] > cmb[depth * 2 - 3] + 1) {	// mathematical shortcut (only valid for depth > 1)
                    cmb[2 * depth - 1] = 2;								        // reset current curl to default
                    depth -= 1;													// and go one back up in depth (we finished current curl/period combination)
                    cmb[depth * 2]++;									        // increase period of this new depth by one
                }
                else if (depth == 1)											// if depth is 1, we can't do anything about 'previous' depths so just break
                    break;
            }

            // in case we changed depths above, we will recalculate the depth to which we will make a combination
            // this uses limits that have been tuned through experiments - there's no theory to this
            if (recalc) {
                int sum = cmb[1] * cmb[1] * cmb[2];                             // weight of first depth
                for (depth = 1; depth < max_depth; depth++) {                   // loop up to max depth
                    if (sum > max_depth * max_depth)
                        break;                                                  // if we go over the limit, we do not use this depth
                    sum += depth * cmb[depth * 2 + 1] * cmb[depth * 2 + 2];   	// otherwise, we continue to next depth
                }
                if (depth == 1 && sum <= ctx.length)                            // counter a couple very difficult edge-case combinations in depth 1
                    depth = 2;
            }
        }
    }

    /*
        once all combinations have been put into the queue, put one last item in the queue
        this has 0 as its first element - this tells the worker to quit once they see this item
    */
    cmb[0] = 0;
    cmbs.push(cmb);

    // update time and log to terminal/file
    now = time(NULL);
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Finished generating combinations\n", &date[4]);
}

/*
    @brief: find the largest power of 2 contained in an (unsigned) int
    input:
    unsigned int n - the number to extract the power from

    returns:
    int power - the power of 2 (not the 2^power value!)
*/
int largestPowerOf2(unsigned int n) {
    int power = 0;
    while (n != 0) {
        n >>= 1;
        power += 1;
    }
    return power;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "Need exactly two CL arguments: length and n_lines" << std::endl;
        return 0;
    }
    
    const int length = atoi(argv[1]);
    const int n_lines = atoi(argv[2]);
    const int max_depth = largestPowerOf2(length);

    char date[26] = {};
    time_t now = time(NULL);
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Started generating combinations for length %i, maximum depth %i\n", &date[4], length, max_depth);

    generate_combinations(length, max_depth, n_lines);

}