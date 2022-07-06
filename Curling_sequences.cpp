// Made by Steven Boonstoppel, with crucial speed improvements thanks to Vladimir Feinstein, algorithm by Levi van de Pol
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 07-06-2021; completed time to length 48: 68 milliseconds*
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, Microsoft VS Studio 2019 Compiler

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
#define FILE std::ofstream file
#define FILE_OPEN file.open("Sequences.txt")
#define OUTPUT stdout // file
#define FILE_CLOSE file.close()
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
std::vector<int> g_best_tails;					// all global longest tail lengths
std::vector<v16_t> g_best_grts;			        // all global generators for the longest tail lengths
std::mutex m_tails, m_queue;					// locks in order to be thread safe when touching the queue and global best_tails/best_grts

// all known record tail lengths so far
std::unordered_map<int, int> known_tails = {
    {2, 2},      {4, 4},      {6, 8},      {8, 58},     {9, 59},     {10, 60},    {11, 112},  {14, 118},
    {19, 119},   {22, 120},   {48, 131},   {68, 132},   {73, 133},   {77, 173},   {85, 179},  {115, 215},
    {116, 228},  {118, 229},  {128, 332},  {132, 340},  {133, 342},  {149, 343},  {154, 356}, {176, 406},
    {197, 1668}, {199, 1669}, {200, 1670}, {208, 1708}, {217, 1836}, {290, 3382}, {385, 3557}
};

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
    @brief: fastest-ish way of erasing an element from a vector:
    find the required value's index through a simple loop
    then swap it with the back of the vector and pop that back value

    input:
    v16_t& vec - the vector from which to remove the value
    int x - the value in question
*/
void erase(v16_t& vec, int x) {
    int i = 0;
    while (vec[i++] != x) {}    // find current index
    vec[i - 1] = vec.back();    // swap with last index
    vec.pop_back();             // remove new last index
}

/*
    @brief: if a curl/period combination has failed, try upgrading the period or else perform backtracking on the tail.
    in the backtracking algorithm, the function `backtracking_step' is the part which actually backtracks.
    this function is called when we come to the conclusion that our current expanded tail does not work.
    An EXPANDED TAIL is defined as the combination of:
    1. the tail values
    2. the MINIMAL corresponding periods.
    therefore, we raise the last element of the expanded tail (which is p_cand) by 1.
    in theory, we could now immediately test this new expanded tail.
    however, to improve efficiency we use a number of tricks to skip options.
    Trick 1. if c_cand * p_cand is larger than the size of the sequence, then this option is impossible.
    Trick 2. For every two consecutive elements a,b of the tail (not the expanded tail) we must have a+1 \geq b
    Trick 3. If the for-last curl has not been changed before, then we dont have to update the for-last period, but we can raise the for-last curl immediately.
    Trick 4. If we raise c_cand, then the corresponding repitition has to go back to the generator of the sequence. This gives a lower bound for p_cand.
    Trick 5. In the case of trick 3, we don't have to change seq and seq_map back.

    Mathematical note: actually, every element of the tail has not just one period, but a set of periods which contains at least one element.
    We don't want to backtrack through all possible sets of periods.
    Therefore, when we say 'period' we mean one of the possible corresponding periods.
    Still, we won't have double cases and we won't miss cases.
    We prove this in our article.
    TODO Explain and prove this in the article. See also the TODO in the function test_seq_new

    input:
    context& ctx - all relevant variables
*/
INLINING
void backtracking_step(context& ctx) {
    ++ctx.p_cand;                                                           // try period one larger now
    while ((ctx.c_cand * ctx.p_cand) > ctx.seq.size()) {                    // trick 1: we look for a pair for which the product is not too large;
                                                                            // as long as this is not the case, we keep backtracking.
        if (ctx.periods.size() < ctx.depth)                                 // stop if tail is smaller than minimum depth; this means that we're completely done.
            break;

        // the product is too high, so c_cand has to be raised (trick 1)
        if (ctx.c_cand <= ctx.seq.back()) {                                 // trick 2 
            ++ctx.c_cand;
            ctx.p_cand = 1 + ctx.periods.size() / ctx.c_cand;              	// trick 4
        }
        else {                                                              // if c_cand > last_curl, then c_cand can not be raised. (trick 2)
                                                                            // therefore, we have to remove the c_cand and p_cand, and change the for-last period.
            int16_t k = (int16_t)ctx.periods.size() - 1;                    // if the for-last curl has never been changed before, we can immediately raise it by 1 (trick 3)
            auto index = ctx.change_indices.find(k);                        // so let's check whether this has been done before
            if (index == ctx.change_indices.end()) {                        // didn't find it?
                ctx.change_indices.insert(k);                               // insert it and let's try increased curl with new period
                ctx.c_cand = ctx.seq.back() + 1;                            // note that ctx.seq.back() will be removed from ctx.seq
                ctx.p_cand = 1 + k / ctx.c_cand;                            // trick 4
                                                                            // trick 5
            }
            else {                                                          // did find it?
                ctx.c_cand = ctx.seq.back();                                // let's try this same curl but now with one higher period
                ctx.p_cand = ctx.periods.back() + 1;
                // we will now change seq and seq_map back to their earlier values. 
                v16_t& temp = ctx.grts_mem[k];                              // retrieve generator from memory. Note that the tail does not have to be 'downgraded', only the generator
                for (int i = 0; i < ctx.length; i++)
                    if (ctx.seq[i] != temp[i]) {                            // check where the differences are, apply them, and change seq_map and seq accordingly
                        erase(ctx.seq_map[ctx.seq[i] + ctx.length], i);
                        ctx.seq_map[temp[i] + ctx.length].push_back(i);
                        ctx.seq[i] = temp[i];
                    }
                ctx.grts_mem.erase(k);                                      // and delete memory entry
            }

            // since we are reducing the length of the tail, remove any entry in the context pointing to the tail element
            auto ind = std::find(ctx.seq_map[ctx.seq.back() + ctx.length].begin(), ctx.seq_map[ctx.seq.back() + ctx.length].end(), ctx.seq.size() - 1);
            ctx.seq_map[ctx.seq.back() + ctx.length].erase(ind);            // remove entry from reverse map
            ctx.seq.pop_back();                                             // and from the sequence
            ctx.periods.pop_back();                                         // and its period
            index = ctx.change_indices.find(k + 1);
            if (index != ctx.change_indices.end())
                ctx.change_indices.erase(index);                            // as well as the index indicating whether we performed backtracking here at one point
        }
    }
}

/*
    @brief: find the true length of the generator for a certain tail.
    a certain tail may not have necessarily required 'length' elements to be set in the generator.
    most often, the first 'n' elements have never been touched for some integer n.
    so for indices i in {0,1,...,n-1}, the value will still be (-length + i).
    those first n elements can then be discarded.

    input:
    context& ctx - all relevant variables

    returns:
    the amount of elements on the right side of the generator that we really need to obtain the tail
*/
int real_generator_length(context& ctx) {
    int i = 0;
    while ((ctx.seq[i] == (-ctx.length + i)) && (++i != ctx.length)) {} // loop through the generator until an index has been changed from its original value
    return (ctx.length - i);                                            // return the total length minus this index, i.e. the true generator length
}

/*
    @brief: in this function we will append the c_cand and p_cand values to the sequence, since they have endured all the tests.

    input:
    context& ctx - all relevant variables
*/
INLINING
void append(context& ctx) {
    for (int i = 0; i < ctx.pairs.size(); i += 2) {             // the sequence changed in test_cands, so it is time to update the map        
        for (int16_t x : ctx.seq_map[ctx.pairs[i + 1] + ctx.length])
            ctx.seq_map[ctx.pairs[i] + ctx.length].push_back(x);// so we move all changed values to their new location
        ctx.seq_map[ctx.pairs[i + 1] + ctx.length].clear();     // and delete their previous entries
                                                                // ctx.pairs is explained in the function test_cands
    }
    ctx.seq.resize(ctx.length);                                 // discard the tail so we can swap ctx.seq into the memory
    ctx.grts_mem[(int16_t)ctx.periods.size()].swap(ctx.seq);
    ctx.seq.swap(ctx.seq_new);                                  // retrieve the new sequence from test_seq_new
    ctx.seq.push_back((int16_t)ctx.c_cand);                     // add the candidates because we passed test_cands and test_seq_new
    ctx.periods.push_back((int16_t)ctx.p_cand);
    int seq_size = (int)ctx.seq.size();
    ctx.seq_map[ctx.c_cand + ctx.length].push_back((int16_t)(seq_size - 1));
    int period = 0;
    while (true) {                                              // as long as the curling number is at least 2,
                                                                // we keep adding the curling number to the sequence.
                                                                // we don't have to perform any tests for this
                                                                // and the sequence doesn't change, so we don't add the index to change_indices
        int curl = krul(ctx.seq, period, seq_size, 2);
        if (curl == 1)
            break;
        ctx.seq.push_back((int16_t)curl);
        ctx.seq_map[curl + ctx.length].push_back((int16_t)seq_size);
        ++seq_size;
        ctx.periods.push_back((int16_t)period);
    }

    int tail = (int)ctx.periods.size();                         // find tail size (seq.size() - length)
    ctx.c_cand = 2;                                             // prepare candidates for the next tests
    ctx.p_cand = 1 + tail / 2;
    ctx.change_indices.insert((int16_t)tail);                   // the sequence will have to change, so we add the tail length to change_indices
    int len = real_generator_length(ctx);                       // retrieve actual generator length
    if (ctx.best_tails[len] < (int16_t)tail) {                  // and update maximum values for this thread
        ctx.best_tails[len] = (int16_t)tail;
        memcpy(&ctx.best_grts[len][0], &ctx.seq[ctx.length - len], sizeof(int16_t) * len);
    }
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

INLINING
void backtracking(context& ctx) {
    if (test_cands(ctx) && test_seq_new(ctx))
        append(ctx);
    else
        backtracking_step(ctx);
}

/*	worker-thread function
    @brief: this function uses backtracking to process all sequences and tails that match the combination it's working on
    it takes a combination from the queue and initiates backtracking to find the largest possible tail for this combination

    input:
    int thread_number - number given to this thread from main for nicer logging
    int len - the user-selected length
*/
void worker(int thread_number, int len) {
    // update time and log to terminal/file
    char date[26] = {};
    time_t now = time(NULL);
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Thread %02d started!\n", &date[4], thread_number);

    // start timer to clock performance
    auto t1 = high_resolution_clock::now();

    /*
        initialize all variables that are required for constructing the curling sequences
        (their explanations are given in the struct)
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

    v16_t cmb;							// initialize empty vector that will hold a combination

    while (true) {
        /*
            lock the queue and copy its exposed front element
            if this element starts with a zero - terminate
            else, pop this element and process it
        */
        {
            std::unique_lock<std::mutex> lock(m_queue);	// lock the queue for thread safety
            if (cmbs.size() == 0)
                continue;								// if the queue is empty, release lock and wait for main thread to fill the queue
            cmb = cmbs.front();				            // if the queue is non-empty, copy its first element
            if (cmb[0] == 0)
                break;									// if the first element in the combination is zero, break and quit
            cmbs.pop();					                // pop combination from the queue
        }

        // store depth as variable in the context
        ctx.depth = cmb[0];

        // intialize an empty sequence with as many elements as the given length
        ctx.seq.resize(ctx.length);

        // fill this sequence as follows: [-n, -n+1, -n+2, ..., -3, -2, -1]
        for (int i = 0; i < ctx.length; ++i)
            ctx.seq[i] = (int16_t)(-ctx.length + i);

        // reset the sequence map
        for (auto& v : ctx.seq_map)
            v.clear();

        // construct the sequence map: this is a reverse of the sequence (value -> index instead of the normal index -> value) for quick accessing
        for (int j = 0; j < ctx.length; ++j)
            ctx.seq_map[j].push_back((int16_t)j);

        // reset the tail (this will hold the periods of all curling numbers in the tail)
        ctx.periods.clear();

        // reset the changed indices (this will hold all indices where we applied backtracking to try and improve)
        ctx.change_indices.clear();

        /*
            simplified version of the curling algorithm that will append 'combination' to the sequence
            extra qualitative explanation can be found in test_cands(), test_seq_new() and append()

            for all curl/period combinations in the current 'combination', change seq, seq_map, periods and change_indices as necessary
        */
        for (int i = 1; i <= ctx.depth; i++) {
            
            ctx.c_cand = cmb[i * 2 - 1];							        // store candidate curl in context
            ctx.p_cand = cmb[i * 2];								        // store candidate period in context

            int period = 0;
            int curl = krul(ctx.seq, period, ctx.length + i - 1, ctx.c_cand);
            if (curl < ctx.c_cand) {
                ctx.change_indices.insert((int16_t)(i - 1));			    // if the curling number is smaller than c_cand,
                                                                            // then the sequence will have to change,
                                                                            // so we add the position to change_indices
                                                                            // note that this is the same thing as what happens in append()
            }
            if (i == ctx.depth)                                             // we will not add the last candidates here,
                break;                                                      // because we are not sure yet whether they can be added

            test_cands(ctx) && test_seq_new(ctx);							// perform test_cands() and test_seq_new() in order to update necessary variables (success by definition)
            for (int j = 0; j < ctx.pairs.size(); j += 2) {
                for (int16_t x : ctx.seq_map[ctx.pairs[j + 1] + ctx.length])
                    ctx.seq_map[ctx.pairs[j] + ctx.length].push_back(x);	// as in append(): update the map with the context that has returned from test_cands() and test_seq_new()
                ctx.seq_map[ctx.pairs[j + 1] + ctx.length].clear();
            }
            ctx.seq.swap(ctx.seq_new);                                                  // as in append(): retrieve the new sequence from test_seq_new
            ctx.seq.push_back((int16_t)ctx.c_cand);                                     // as in append(): add the curl
            ctx.periods.push_back((int16_t)ctx.p_cand);						            // as in append(): add the period
            ctx.seq_map[ctx.c_cand + ctx.length].push_back((int16_t)(ctx.length + 1));	// as in append(): add the curl to the reverse map
        }

        /*
            finally, we will start backtracking using the last depth of 'combination'
            as long as the length of the tail is larger than the used depth, we will set another step in the algorithm
            this will either result in a longer tail, or a backtracking step on the tail and making a small change
            rinse and repeat
        */
     
        do backtracking(ctx);
        while (ctx.periods.size() >= ctx.depth);                            // perform backtracking for this combination (c_cand, p_cand)
    }

    /*
        once the worker has found a zero in the queue, it exits the loop
        and synchronizes its best found tails to the global vectors under a thread safe lock
    */
    std::unique_lock<std::mutex> lock(m_tails);                             // lock global vectors
    for (int i = 0; i <= ctx.length; ++i) {
        if (ctx.best_tails[i] > g_best_tails[i]) {							// if this thread has found a record tail, update value and copy its generator
            g_best_tails[i] = ctx.best_tails[i];
            memcpy(&g_best_grts[i][0], &ctx.best_grts[i][0], sizeof(int16_t) * ctx.length);
        }
    }

    // update time and log to terminal/file
    now = time(NULL);
    UPDATE_TIME;
    auto t2 = high_resolution_clock::now();
    int elapsed = (int)duration_cast<milliseconds>(t2 - t1).count();
    fprintf(OUTPUT, "[%.15s] Thread %02d finished, duration: %ims\n", &date[4], thread_number, elapsed);
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
void generate_combinations(int len, int max_depth) {
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
        while (cmbs.size() < ctx.length * ctx.length && cmb[1] <= ctx.length) {
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
                std::unique_lock<std::mutex> lock(m_queue);
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

        /*
            sleep a little to prevent clogging the queue lock and use too much resources
        */
        std::this_thread::sleep_for(milliseconds(100));
    }

    /*
        once all combinations have been put into the queue, put one last item in the queue
        this has 0 as its first element - this tells the worker to quit once they see this item
    */
    cmb[0] = 0;
    cmbs.push(cmb);

    // update time and log to terminal/file
    char date[26] = {};
    time_t now = time(NULL);
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Finished generating combinations\n", &date[4]);
}

/*
    @brief: log the results
    and check them against known record values to quickly see if a difference has been found

    input:
    int max_depth - the calculated maximum depth for the user-defined length
*/
void log_results(int max_depth) {

    //update time and log to terminal/file
    time_t now = time(NULL);
    char date[26] = {};
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Logging results:\n", &date[4]);

    /*
        find all positions where the maximum tail length is increased
        these are the record strings of interest
        for ease of use, compare them to known record values and write to terminal/file
    */
    int record = 0;
    for (int i = 0; i < g_best_tails.size(); ++i)
        if (g_best_tails[i] > record) {
            record = g_best_tails[i];
            if (i <= max_depth)
                continue;
            if (known_tails.find(i) == known_tails.end())
                fprintf(OUTPUT, "NEW: ");
            else if (known_tails[i] != record)
                fprintf(OUTPUT, "???: ");
            else
                fprintf(OUTPUT, "OLD: ");
            fprintf(OUTPUT, "%i: %i, [", i, record);
            for (int x : g_best_grts[i])
                if (x != 0)
                    fprintf(OUTPUT, "%i,", x);
            fprintf(OUTPUT, "\b]\n");
        }
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
        std::cout << "Need exactly two CL arguments: length and number of threads" << std::endl;
        //return 0;
    }
    
    const int length = 100;// atoi(argv[1]);
    const int thread_count = 16;// atoi(argv[2]);
    const int max_depth = largestPowerOf2(length);

    char date[26] = {};
    time_t now = time(NULL);
    UPDATE_TIME;
    fprintf(OUTPUT, "[%.15s] Started length %i, maximum depth %i, thread count %i\n", &date[4], length, max_depth, thread_count);

    g_best_tails.resize(length + 1);
    g_best_grts.resize(length + 1);
    for (int i = 0; i <= length; i++)
        g_best_grts[i].resize(length);

    FILE_OPEN;

    // Start the calculations
    std::vector<std::thread> thread_vector;
    thread_vector.emplace_back(std::thread(generate_combinations, length, max_depth));
    for (int thread_number = 1; thread_number <= thread_count; ++thread_number)
        thread_vector.emplace_back(std::thread(worker, thread_number, length));
    for (auto& th : thread_vector)
        th.join();

    // Log the results
    log_results(max_depth);

    FILE_CLOSE;
}