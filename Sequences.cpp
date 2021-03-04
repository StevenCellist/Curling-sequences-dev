// Made by Steven Boonstoppel, with crucial speed improvements thanks to StackOverflow user Vlad Feinstein
// First version: 18-11-2020; estimated time to length 48: 250 years*
// Current version: 11-02-2021; completed time to length 48: 3 seconds*
// * reference CPU: AMD Ryzen 7 3800X, 16 threads @ ~4.2 GHz boost, Microsoft VS Studio 2019 Compiler

#include <vector>
#include <array>
#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>
#include <fstream>
#include <algorithm>
#include <immintrin.h>
#include <set>
using namespace std::chrono;

std::mutex m_tails;						// lock to avoid conflicting changes while editing from multiple threads
const int thread_count = std::thread::hardware_concurrency();

const int length = 90;					// length to which the program will run
const int fraction = 4;					// FRACTION OF THE SEQUENCE (MATHEMATICALLY INCORRECT, AFFECTS RESULTS FOR > 4)
const int maximum = 254;				// MAXIMUM NUMBER OF BITS BEFORE TRAP (I got BSOD when I tried 255 once, IDK if that was due to the program)
const int cache_bits = 30;              // cache is built up to this length
int global_tails[length + 1];			// array to hold the maximum tail lengths for each length
thread_local int tails[length + 1] = {};// array to hold the maximum tail lengths for each length, per thread

uint8_t* cache[cache_bits + 1];         // allocate starting with cache[0], although it's not used, to make indexing faster
uint8_t* cache30;						// fast pointer to the cache position 'cache_bits'

#pragma region SSE/AVX2 LSFRs and compares, small cache
struct m512 {
	m512(__m256i lo = _mm256_setzero_si256(), __m256i hi = _mm256_setzero_si256()) : m_lo(lo), m_hi(hi) {}
	__m256i m_lo, m_hi;
};

// main sources: 
// http://notabs.org/lfsr/software/
// https://stackoverflow.com/questions/9980801/looking-for-sse-128-bit-shift-operation-for-non-immediate-shift-value

// bit shift right a 128-bit value
// __m128i data - data to shift
//    int count - number of bits to shift (< 64)
__m128i shift_right(__m128i data, const int count) {
	__m128i carry = _mm_bsrli_si128(data, 8);       // Shift the whole thing 64 bits right
	if (count > 63) {
		return _mm_srli_epi64(carry, count - 64);
	}
	carry = _mm_slli_epi64(carry, 64 - count);		// After bslri shifted right by 64b
	data = _mm_srli_epi64(data, count);
	return _mm_or_si128(data, carry);
}

// bit shift left a 128-bit value by one bit
// __m128i* data - data to shift
void shift_left_one(__m128i* data) {
	__m128i carry = _mm_bslli_si128(*data, 8);      // Shift the whole thing 64 bits right
	carry = _mm_srli_epi64(carry, 63);				// After bslli shifted left by 64b
	*data = _mm_slli_epi64(*data, 1);
	*data = _mm_or_si128(*data, carry);
}

void shift_left(__m128i* data, const int count) {
	__m128i carry = _mm_bslli_si128(*data, 8);      // Shift the whole thing 64 bits right
	if (count > 63) {
		*data = _mm_slli_epi64(carry, count - 64);
	}
	else {
		carry = _mm_srli_epi64(carry, 64 - count);		// After bslli shifted left by 64b
		*data = _mm_slli_epi64(*data, count);
		*data = _mm_or_si128(*data, carry);
	}
}

void shift_right(__m256i* data, int count) {
	if (count >= 192) {
		*data = _mm256_permute4x64_epi64(*data, 0x1B);						// order: 00 01 10 11
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0x03);	// clear 3 upper qwords
		count -= 192;
	}
	else if (count >= 128) {
		*data = _mm256_permute4x64_epi64(*data, 0x1E);						// order: 00 01 11 10
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0x0F);	// clear 2 upper qwords
		count -= 128;
	}
	else if (count >= 64) {
		*data = _mm256_permute4x64_epi64(*data, 0x39);						// order: 00 11 10 01
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0x3F);	// clear 1 upper qword
		count -= 64;
	}
	if (count) {
		__m256i innerCarry = _mm256_slli_epi64(*data, 64 - count);				// carry outs in bit 0 of each qword
		__m256i rotate = _mm256_permute4x64_epi64(innerCarry, 0x39);			// rotate ymm left 64 bits
		innerCarry = _mm256_blend_epi32(_mm256_setzero_si256(), rotate, 0x3F);	// clear upper qword
		*data = _mm256_srli_epi64(*data, count);								// shift all qwords right
		*data = _mm256_or_si256(*data, innerCarry);								// propagate carrys from low qwords
	}
}

// bit shift left by one a 256-bit value using ymm registers
//          __m256i* data - data to shift
void shift_left_one(__m256i* data) {
	__m256i innerCarry = _mm256_srli_epi64(*data, 63);                     // carry outs in bit 0 of each qword
	__m256i rotate = _mm256_permute4x64_epi64(innerCarry, 0x93);           // rotate ymm left 64 bits
	innerCarry = _mm256_blend_epi32(_mm256_setzero_si256(), rotate, 0xFC); // clear lower qword
	*data = _mm256_slli_epi64(*data, 1);                                   // shift all qwords left
	*data = _mm256_or_si256(*data, innerCarry);                            // propagate carrys from low qwords
}
// bit shift left a 256-bit value using ymm registers
//          __m256i* data - data to shift
void shift_left(__m256i* data, int count) {
	if (count >= 192) {														// only 1 left qword matter
		*data = _mm256_permute4x64_epi64(*data, 0x1B);						// order: 00 01 10 11
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0xC0);	// clear 3 lower qwords
		count -= 192;
	}
	else if (count >= 128) {												// only 2 left qwords matter
		*data = _mm256_permute4x64_epi64(*data, 0x4E);						// order: 01 00 11 10
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0xF0);	// clear 2 lower qwords
		count -= 128;
	}
	else if (count >= 64) {													// 3 left qwords matter
		*data = _mm256_permute4x64_epi64(*data, 0x93);						// order: 10 01 00 11
		*data = _mm256_blend_epi32(_mm256_setzero_si256(), *data, 0xFC);	// clear 1 lower qword
		count -= 64;
	}
	if (count) {
		__m256i innerCarry = _mm256_srli_epi64(*data, 64 - count);             // carry outs in bit 0 of each qword
		__m256i rotate = _mm256_permute4x64_epi64(innerCarry, 0x93);           // rotate ymm left 64 bits
		innerCarry = _mm256_blend_epi32(_mm256_setzero_si256(), rotate, 0xFC); // clear lower qword
		*data = _mm256_slli_epi64(*data, count);                               // shift all qwords left
		*data = _mm256_or_si256(*data, innerCarry);                            // propagate carrys from low qwords
	}
}

// bit shift left by one a 256-bit value using ymm registers
//          __m256i* data - data to shift
void push(m512* data, bool one) {
	static const __m256i mask1 = _mm256_set_epi64x(0, 0, 0, 1);
	__m256i carry = data->m_lo;
	shift_right(&carry, 256 - 1);						// get carry bits
	shift_left_one(&data->m_lo);						// shift low part
	shift_left_one(&data->m_hi);						// shift high part
	data->m_hi = _mm256_or_si256(data->m_hi, carry);	// propagate carrys from m_lo to m_hi
	if (one)
		data->m_lo = _mm256_or_si256(data->m_lo, mask1);// and add 1 if the curl is 3 (binary twin)
}
// bit shift left by one a 256-bit value using ymm registers
//          __m256i* data - data to shift
void shift_left(m512* data, int count) {
	__m256i carry = data->m_lo;
	shift_right(&carry, 256 - count);					// get carry bits
	shift_left(&data->m_lo, count);						// shift low part
	shift_left(&data->m_hi, count);						// shift high part
	data->m_hi = _mm256_or_si256(data->m_hi, carry);	// propagate carrys from m_lo to m_hi
}

void shift_right(m512* data, int count) {
	__m256i carry = data->m_hi;
	shift_left(&carry, 256 - count);					// get carry bits
	shift_right(&data->m_lo, count);					// shift low part
	shift_right(&data->m_hi, count);					// shift high part
	data->m_lo = _mm256_or_si256(data->m_lo, carry);	// propagate carrys from m_lo to m_hi
}

// check if 2 parts of the sequence are equal
//    __m128i* src - last part of the sequence to compare
// __m128i* target - shifted part of the sequence to compare
//   __m128i* mask - mask for the part of the sequence we are interested in
bool masked_eq(const __m128i* src, const __m128i* target, const __m128i* mask) {
	__m128i diff = _mm_xor_si128(*src, *target);
	return (_mm_test_all_zeros(diff, *mask));
}

// check if 2 parts of the sequence are equal
//    __m256i* src - last part of the sequence to compare
// __m256i* target - shifted part of the sequence to compare
//   __m256i* mask - mask for the part of the sequence we are interested in
bool masked_eq(const __m256i* src, const __m256i* target, const __m256i* mask) {
	__m256i diff = _mm256_xor_si256(*src, *target);
	return (_mm256_testz_si256(diff, *mask));
}

// mask of 64 bits to select last bits of a sequence
const uint64_t mask64[65] = { 0,
  0x1, 0x3, 0x7, 0xF, 0x1F, 0x3F, 0x7F, 0xFF,
  0x1FF, 0x3FF, 0x7FF, 0xFFF, 0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF,
  0x1FFFF, 0x3FFFF, 0x7FFFF, 0xFFFFF, 0x1FFFFF, 0x3FFFFF, 0x7FFFFF, 0xFFFFFF,
  0x1FFFFFF, 0x3FFFFFF, 0x7FFFFFF, 0xFFFFFFF, 0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF,
  0x1FFFFFFFF, 0x3FFFFFFFF, 0x7FFFFFFFF, 0xFFFFFFFFF, 0x1FFFFFFFFF, 0x3FFFFFFFFF, 0x7FFFFFFFFF, 0xFFFFFFFFFF,
  0x1FFFFFFFFFF, 0x3FFFFFFFFFF, 0x7FFFFFFFFFF, 0xFFFFFFFFFFF, 0x1FFFFFFFFFFF, 0x3FFFFFFFFFFF, 0x7FFFFFFFFFFF, 0xFFFFFFFFFFFF,
  0x1FFFFFFFFFFFF, 0x3FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0xFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFF, 0x3FFFFFFFFFFFFF, 0x7FFFFFFFFFFFFF, 0xFFFFFFFFFFFFFF,
  0x1FFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF,
};
const uint64_t cache_mask = mask64[cache_bits];

// mask of 128 bits to select last bits of a sequence
const std::array<__m128i, 129> build_mask() {
	std::array<__m128i, 129> mask128;
	for (int i = 0; i < 128; ++i)
		mask128[i] = i < 64 ?
		_mm_set_epi64x(0, mask64[i % 64]) :
		_mm_set_epi64x(mask64[i % 64], mask64[64]);

	mask128[128] = _mm_set_epi64x(mask64[64], mask64[64]);
	return mask128;
}
const std::array<__m128i, 129> mask128 = build_mask();

// mask of 256 bits to select last bits of a sequence
const std::array<__m256i, 257> build_mask_256() {
	std::array<__m256i, 257> mask;
	for (int i = 0; i < 256; ++i)
		mask[i] = i < 128 ?
		_mm256_set_m128i(mask128[0], mask128[i % 128]) :
		_mm256_set_m128i(mask128[i % 128], mask128[128]);

	mask[256] = _mm256_set_m128i(mask128[128], mask128[128]);
	return mask;
}
const std::array<__m256i, 257> mask256 = build_mask_256();

// pre-calculate the upper bounds (for pattern_length) of finding certain curls for all lengths
const std::array< std::array<uint16_t, 4>, 512> build_limit_cache() {
	std::array< std::array<uint16_t, 4>, 512> limit_cache;
	for (uint16_t i = 0; i < 512; ++i) {
		for (uint8_t c = 1; c < 4; ++c) { // ignore index 1 for easy access by curl 1..3
			limit_cache[i][c] = i / (c + 1);
		}
	}
	return limit_cache;
}
const std::array< std::array<uint16_t, 4>, 512> limit_cache = build_limit_cache();

// pre-calculate the lower bound on pattern lengths
const uint8_t pattern_length_cache[4] = { 0, cache_bits / 2 + 1, cache_bits / 3 + 1, cache_bits / 4 + 1 };

struct OptLen {
	__m256i opt;
	uint8_t len;
	OptLen(__m256i opt, uint8_t len) : opt(opt), len(len) {}
};
struct compare_OptLen {
	bool operator()(const OptLen& a, const OptLen& b) const {
		if (a.len < b.len)
			return true;
		if (a.len > b.len)
			return false;
		for (int i = 3; i >= 0; --i) {
			if (a.opt.m256i_u64[i] < b.opt.m256i_u64[i])
				return true;
			if (a.opt.m256i_u64[i] > b.opt.m256i_u64[i])
				return false;
		}
		return false;
	}
};

#pragma endregion

#pragma region SEQUENCE BUILDING
// krul function to calculate curling number for sequence length <= 64
//		 const uint8_t i - length of sequence
//   bool last_curl_is_2 - value of previous curling number (binary twin)
//  const uint64_t seq64 - current sequence (used for <= 64 bits)
__forceinline
uint8_t krul64(const uint8_t i, bool last_curl_is_2, const uint64_t seq64) {

	// extract curl from cache
	uint64_t key = seq64 & cache_mask;								// find the 'key' of the sequence for the cached amount of bits
	uint8_t curl = cache30[key];									// retrieve the curl that matches this key
	key = (key * 2) & cache_mask;									// the next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch next needed part of cache from RAM to L1 cache
	key = (key * 2) & cache_mask;									// the second next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch second next needed part of cache from RAM to L1 cache

	if (curl > (3 - last_curl_is_2))								// if the last curl was 2 resp. 3 and the current is 3 resp. 4, we can break (mathematical)
		return curl;

	uint8_t pattern_length = pattern_length_cache[curl];			// pattern_length is the length of the period we check, starting from the lower bound
	const uint16_t limit_4 = limit_cache[i][3];						// we can find a curl of 4 only up to limit_4
	const uint16_t limit_3 = limit_cache[i][2];						// we can find a curl of 3 only up to limit_3
	const uint16_t limit_2 = limit_cache[i][1];						// we can find a curl of 2 only up to limit_2

	// check the sequence for repetiton of up to 4 times
	for (; pattern_length <= limit_4; ++pattern_length) {
		uint64_t mask = mask64[pattern_length];						// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			target >>= pattern_length;								// shift the target sequence by pattern_length bits
			if (!((target ^ seq64) & mask)) {						// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				target >>= pattern_length;							// shift the copy of sequence by pattern_length bits
				if (!((target ^ seq64) & mask))						// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
			else {
				if (curl == 1)
					curl = 2;
			}
		}
	}
	if (curl == 3)
		return curl;												// we won't find curl > 3 from now on because we passed a period length of 1/4th of the sequence

	// check the sequence for repetition of up to 3 times
	for (; pattern_length <= limit_3; ++pattern_length) {
		uint64_t mask = mask64[pattern_length];						// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask)) {							// continue if the parts match
			target >>= pattern_length;								// shift the target sequence by pattern_length bits
			if (!((target ^ seq64) & mask))							// check if they match; if yes, finish, because we won't find frequency > 3
				return 3;
			if (curl == 1) {
				curl = 2;
			}
		}
	}
	if (curl == 2)
		return curl;												// we won't find curl > 2 from now on because we passed a period length of 1/3th of the sequence

	// check the sequence for repetition of up to 3 times
	for (; pattern_length <= limit_2; ++pattern_length) {
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask64[pattern_length]))			// check if they match; if yes, finish, because we won't find frequency > 2
			return 2;
	}
	return curl;													// no repetition found
}

// krul function to calculate curling number for 64 <= sequence length <= 128
//		 const uint8_t i - length of sequence
//   bool last_curl_is_2 - value of previous curling number (binary twin)
// const __m128i* seq128 - current sequence (used for <= 128 bits)
__forceinline
uint8_t krul128(const uint8_t i, bool last_curl_is_2, const __m128i* seq128) {

	// extract curl from cache
	uint64_t key = seq128->m128i_u64[0] & cache_mask;				// find the 'key' of the sequence for the cached amount of bits
	uint8_t curl = cache30[key];									// retrieve the curl that matches this key
	key = (key * 2) & cache_mask;									// the next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch next needed part of cache from RAM to L1 cache
	key = (key * 2) & cache_mask;									// the second next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch second next needed part of cache from RAM to L1 cache

	if (curl > (3 - last_curl_is_2))								// if the last curl was 2 resp. 3 and the current is 3 resp. 4, we can break (mathematical)
		return curl;

	uint8_t pattern_length = pattern_length_cache[curl];			// pattern_length is the length of the period we check, starting from the lower bound
	const uint16_t limit_3 = limit_cache[i][2];						// we can find a curl of 3 only up to limit_3
	const uint16_t limit_2 = limit_cache[i][1];						// we can find a curl of 2 only up to limit_2

	const uint64_t seq64 = seq128->m128i_u64[0];					// copy of last 64 bits of the sequence for faster curl calculation

	// check the sequence for repetition of up to 4 times
	for (; pattern_length <= 16; ++pattern_length) {				// optimized loop for pattern_length <= 16
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			target >>= pattern_length;								// shift the target sequence by pattern_length bits
			if (!((target ^ seq64) & mask)) {						// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				target >>= pattern_length;							// shift the copy of sequence by pattern_length bits
				if (!((target ^ seq64) & mask))						// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
			else {
				if (curl == 1)
					curl = 2;
			}
		}
	}

	// check the sequence for repetition of up to 4 times or up to limit
	int current_limit = std::min(limit_cache[i][curl], (uint16_t)32);
	for (; pattern_length <= current_limit; ++pattern_length) {		// optimized loop for pattern_length <= 32, or we stop earlier if limit < 32
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits (short shift)
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			if (curl == 1) {
				curl = 2;
				int limit = limit_cache[i][curl];					// update the limit if a new curl is found
				if (current_limit > limit)
					current_limit = limit;							// update the current_limit if the new limit is lower
			}
			__m128i t128 = shift_right(*seq128, pattern_length * 2);	// shift the target sequence twice by pattern_length bits (long shift because we cross 64 bits)
			if (!((t128.m128i_u64[0] ^ seq64) & mask)) {			// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				t128 = shift_right(t128, pattern_length);		// shift the copy of sequence by pattern_length bits
				if (!((t128.m128i_u64[0] ^ seq64) & mask))			// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				if (curl == 2) {
					curl = 3;
					int limit = limit_cache[i][curl];				// update the limit if a new curl is found
					if (current_limit > limit)
						current_limit = limit;						// update the current_limit if the new limit is lower
				}
			}
		}
	}
	if (curl == 3)
		return curl;												// we won't find curl > 3 from now on because we passed a period length of 1/4th of the sequence

	// check the sequence for repetition of up to 3 times
	for (; pattern_length <= limit_3; ++pattern_length) {			// optimized loop for pattern_length <= 64
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t pattern = seq128->m128i_u64[0] & mask;				// copy of original sequence
		__m128i temp = shift_right(*seq128, pattern_length);  // shift the copy of sequence by pattern_length bits
		if ((temp.m128i_u64[0] & mask) == pattern) {				// check if they match; if yes, continue (frequency = 2)
			if (curl == 1) {
				curl = 2;
			}
			temp = shift_right(temp, pattern_length);         // shift the copy of sequence by pattern_length bits
			if ((temp.m128i_u64[0] & mask) == pattern) {			// check if they match; if yes, finish, because we won't find frequency > 3
				return 3;
			}
		}
	}
	if (curl == 2)                                                  // we won't find curl > 2 from now on because we passed a period length of 1/3th of the sequence
		return curl;

	// check the sequence for repetition of up to 2 times
	for (; pattern_length <= limit_2; ++pattern_length) {			// optimized loop for pattern_length <= 64
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t pattern = seq128->m128i_u64[0] & mask;				// masked copy of original sequence
		__m128i temp = shift_right(*seq128, pattern_length);  // shift the copy of sequence by pattern_length bits
		if ((temp.m128i_u64[0] & mask) == pattern) {				// check if they match; if yes, finish, because we won't find frequency > 2
			return 2;
		}
	}

	return curl;													// no repetition found
}

// krul function to calculate curling number for 128 <= sequence length <= 250
//		 const uint8_t i - length of sequence
//   bool last_curl_is_2 - value of previous curling number (binary twin)
// const __m256i* seq256 - current sequence (used for <= 256 bits)
__forceinline
uint8_t krul256(const uint8_t i, bool last_curl_is_2, const __m256i* seq256) {

	// extract curl from cache
	uint64_t key = seq256->m256i_u64[0] & cache_mask;				// find the 'key' of the sequence for the cached amount of bits
	uint8_t curl = cache30[key];									// retrieve the curl that matches this key
	key = (key * 2) & cache_mask;									// the next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch next needed part of cache from RAM to L1 cache
	key = (key * 2) & cache_mask;									// the second next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch second next needed part of cache from RAM to L1 cache

	if (curl > (3 - last_curl_is_2))								// if the last curl was 2 resp. 3 and the current is 3 resp. 4, we can break (mathematical)
		return curl;

	uint8_t pattern_length = pattern_length_cache[curl];			// pattern_length is the length of the period we check, starting from the lower bound

	__m128i seq128 = _mm256_castsi256_si128(*seq256);				// cast the low half of 256-bit sequence to 128-bit for faster functions
	uint64_t seq64 = seq256->m256i_u64[0];							// cast the low quarter of 256-bit sequence to 64-bit for faster functions

	// check the sequence for repetition of up to 4 times
	for (; pattern_length <= 16; ++pattern_length) {				// optimized loop for pattern_length <= 16
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			target >>= pattern_length;								// shift the copy of sequence by pattern_length bits
			if (!((target ^ seq64) & mask)) {						// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				target >>= pattern_length;							// shift the copy of sequence by pattern_length bits
				if (!((target ^ seq64) & mask))						// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
			else {
				if (curl == 1)
					curl = 2;
			}
		}
	}

	// check the sequence for repetition of up to 4 times
	for (; pattern_length <= 32; ++pattern_length) {				// optimized loop for pattern_length <= 32
		__m256i temp = *seq256;										// copy of original sequence
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits (short shift)
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			if (curl == 1)
				curl = 2;
			__m128i t128 = shift_right(_mm256_castsi256_si128(temp), pattern_length * 2);	// shift the target sequence twice by pattern_length bits (because we cross 64 bits)
			if (!((t128.m128i_u64[0] ^ seq64) & mask)) {			// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				t128 = shift_right(t128, pattern_length);		// shift the copy of sequence by pattern_length bits
				if (!((t128.m128i_u64[0] ^ seq64) & mask))			// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
		}
	}

	// check the sequence for repetition of up to 4 times or up to limit
	int limit = limit_cache[i][curl];
	int current_limit = std::min(limit, 63);
	for (; pattern_length <= current_limit; ++pattern_length) {		// optimized loop for pattern_length <= 63
		__m256i temp = *seq256;										// copy of original sequence
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t pattern = seq256->m256i_u64[0] & mask;				// masked copy of original sequence
		shift_right(&temp, pattern_length);							// shift the copy of sequence by pattern_length bits
		if ((temp.m256i_u64[0] & mask) == pattern) {				// check if they match; if yes, continue (frequency = 2)
			if (curl == 1) {
				curl = 2;
				limit = limit_cache[i][curl];						// update the limit if a new curl is found
				if (current_limit > limit)
					current_limit = limit;							// update the current_limit if the new limit is lower
			}
			int matching_length = 3 * pattern_length;
			if (matching_length <= i) {								// check if another period fits in the sequence length
				shift_right(&temp, pattern_length);					// shift the copy of sequence by pattern_length bits
				if ((temp.m256i_u64[0] & mask) == pattern) {		// check if they match; if yes, continue (frequency = 2)
					if (last_curl_is_2)								// if the last curl was 2, we can break (mathematical)
						return 3;
					matching_length += pattern_length;
					if (matching_length <= i) {						// check if another period fits in the sequence length
						shift_right(&temp, pattern_length);			// shift the copy of sequence by pattern_length bits
						if ((temp.m256i_u64[0] & mask) == pattern) {// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
							return 4;
						}
					}
					if (curl == 2) {
						curl = 3;
						limit = limit_cache[i][curl];				// update the limit if a new curl is found
						if (current_limit > limit)
							current_limit = limit;					// update the current_limit if the new limit is lower
					}
				}
			}

		}
	}

	if (limit < pattern_length)
		return curl;

	if (curl == 3)
		return curl;												// we won't find curl > 3 from now on because we passed a period length of 1/4th of 250

	// optimized case for 64 bit pattern_length

	if (seq256->m256i_u64[0] == seq256->m256i_u64[1]) {				// check if the next 64-bit part matches; if yes, continue (frequency = 2)
		if (i >= 192 && seq256->m256i_u64[0] == seq256->m256i_u64[2]) {	// check if the next 64-bit part matches; if yes, finish, because we won't find frequency > 3
			return 3;
		}
		if (curl == 1) {
			curl = 2;
			limit = limit_cache[i][curl];							// update the limit if a new curl is found
		}
	}
	pattern_length++;

	current_limit = std::min(limit, 85);
	for (; pattern_length <= current_limit; ++pattern_length) {     // optimized loop for pattern_length <= 128
		__m256i temp = *seq256;										// copy of original sequence
		const __m128i mask = mask128[pattern_length];				// mask for this pattern_length
		shift_right(&temp, pattern_length);						// shift the copy of sequence by pattern_length bits
		__m128i t128 = _mm256_castsi256_si128(temp);				// cast the low half of 256-bit sequence to 128-bit for faster functions
		if (masked_eq(&t128, &seq128, &mask)) {                   // check if they match; if yes, continue (frequency = 2)
			shift_right(&temp, pattern_length);					// shift the copy of sequence by pattern_length bits
			t128 = _mm256_castsi256_si128(temp);					// cast the low half of 256-bit sequence to 128-bit for faster functions
			if (masked_eq(&t128, &seq128, &mask))                 // check if they match; if yes, finish, because we won't find frequency > 3
				return 3;
			if (2 > curl) {
				curl = 2;
				limit = limit_cache[i][curl];						// update the limit if a new curl is found
				if (current_limit > limit)
					current_limit = limit;							// update the current_limit if the new limit is lower
			}
		}
	}
	if (curl == 2)                                                  // we won't find curl > 2 from now on because we passed a period length of 1/3th of 250
		return curl;

	for (; pattern_length <= limit; ++pattern_length) {             // optimized loop for pattern_length <= 128
		__m256i temp = *seq256;										// copy of original sequence
		shift_right(&temp, pattern_length);						// shift the copy of sequence by pattern_length bits
		__m128i t128 = _mm256_castsi256_si128(temp);				// cast the low half of 256-bit sequence to 128-bit for faster functions
		if (masked_eq(&t128, &seq128, &mask128[pattern_length]))  // check if they match; if yes, finish, because we won't find frequency > 2
			return 2;
	}
	return curl;													// no repetition found
}

// krul function to calculate curling number for 255 <= sequence length <= 512
//		 const uint8_t i - length of sequence
//   bool last_curl_is_2 - value of previous curling number (binary twin)
//    const m512* seq512 - current sequence (used for <= 512 bits)
//__forceinline
uint8_t krul512(const uint16_t i, bool last_curl_is_2, const m512 const* seq512) {

	// extract curl from cache
	uint64_t key = seq512->m_lo.m256i_u64[0] & cache_mask;			// find the 'key' of the sequence for the cached amount of bits
	uint8_t curl = cache30[key];									// retrieve the curl that matches this key
	key = (key * 2) & cache_mask;									// the next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch next needed part of cache from RAM to L1 cache
	key = (key * 2) & cache_mask;									// the second next sequence will likely be shifted left
	_mm_prefetch((const char*)&cache30[key], _MM_HINT_T0);			// prefetch second next needed part of cache from RAM to L1 cache

	if (curl > (3 - last_curl_is_2)) 								// if the last curl was 2 resp. 3 and the current is 3 resp. 4, we can break (mathematical)
		return curl;

	uint16_t pattern_length = pattern_length_cache[curl];			// pattern_length is the length of the period we check, starting from the lower bound

	__m128i seq128 = _mm256_castsi256_si128(seq512->m_lo);			// cast the low half of 256-bit sequence to 128-bit for faster functions
	uint64_t seq64 = seq512->m_lo.m256i_u64[0];						// cast the low quarter of 256-bit sequence to 64-bit for faster functions

	// check the sequence for repetition of up to 4 times
	for (; pattern_length <= 16; ++pattern_length) {				// optimized loop for pattern_length <= 16
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			target >>= pattern_length;								// shift the copy of sequence by pattern_length bits
			if (!((target ^ seq64) & mask)) {						// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				target >>= pattern_length;							// shift the copy of sequence by pattern_length bits
				if (!((target ^ seq64) & mask)) 					// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
			else {
				if (curl == 1)
					curl = 2;
			}
		}
	}

	// check the sequence for repetition of up to 4 times
	for (; pattern_length <= 32; ++pattern_length) {				// optimized loop for pattern_length <= 32
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t target = seq64 >> pattern_length;					// copy of original sequence, shifted by pattern_length bits (short shift)
		if (!((target ^ seq64) & mask)) {							// check if they match; if yes, continue (frequency = 2)
			if (curl == 1)
				curl = 2;
			__m128i t128 = shift_right(seq128, pattern_length * 2);	// shift the target sequence twice by pattern_length bits (because we cross 64 bits)
			if (!((t128.m128i_u64[0] ^ seq64) & mask)) {			// check if they match; if yes, continue (frequency = 3)
				if (last_curl_is_2) 								// if the last curl was 2, we can break (mathematical)
					return 3;
				t128 = shift_right(t128, pattern_length);			// shift the copy of sequence by pattern_length bits
				if (!((t128.m128i_u64[0] ^ seq64) & mask))			// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				curl = 3;
			}
		}
	}

	for (; pattern_length <= 63; ++pattern_length) {				// optimized loop for pattern_length <= 63
		__m256i temp = seq512->m_lo;								// copy of original sequence
		const uint64_t mask = mask64[pattern_length];				// mask for this pattern_length
		uint64_t pattern = seq512->m_lo.m256i_u64[0] & mask;		// masked copy of original sequence
		shift_right(&temp, pattern_length);							// shift the copy of sequence by pattern_length bits
		if ((temp.m256i_u64[0] & mask) == pattern) {				// check if they match; if yes, continue (frequency = 2)
			if (curl == 1)
				curl = 2;
			shift_right(&temp, pattern_length);						// shift the copy of sequence by pattern_length bits
			if ((temp.m256i_u64[0] & mask) == pattern) {			// check if they match; if yes, continue (frequency = 2)
				if (last_curl_is_2)									// if the last curl was 2, we can break (mathematical)
					return 3;
				shift_right(&temp, pattern_length);					// shift the copy of sequence by pattern_length bits
				if ((temp.m256i_u64[0] & mask) == pattern) {		// check if they match; if yes, finish, because frequency = 4 already terminates the sequence
					return 4;
				}
				if (curl == 2)
					curl = 3;
			}
		}
	}

	// optimized case for 64 bit pattern_length
	if (seq512->m_lo.m256i_u64[0] == seq512->m_lo.m256i_u64[1]) {	// check if the next 64-bit part matches; if yes, continue (frequency = 2)
		if (i >= 192 &&
			seq512->m_lo.m256i_u64[0] == seq512->m_lo.m256i_u64[2]) {	// check if the next 64-bit part matches;
			if (i >= 256 &&
				seq512->m_lo.m256i_u64[0] == seq512->m_lo.m256i_u64[3]) {	// check if the next 64-bit part matches;
				return 4;
			}
			curl = 3;
		}
		if (curl == 1) {
			curl = 2;
		}
	}
	pattern_length++;

	for (; pattern_length <= limit_cache[i][curl]; ++pattern_length) { // optimized loop for pattern_length <= 128
		m512 temp = *seq512;										// copy of original sequence
		const __m256i mask = mask256[pattern_length];				// mask for this pattern_length
		shift_right(&temp, pattern_length);							// shift the copy of sequence by pattern_length bits
		if (masked_eq(&temp.m_lo, &seq512->m_lo, &mask)) {          // check if they match; if yes, continue (frequency = 2)
			if (2 > curl) {
				curl = 2;
			}
			shift_right(&temp, pattern_length);						// shift the copy of sequence by pattern_length bits
			if (masked_eq(&temp.m_lo, &seq512->m_lo, &mask)) {      // check if they match; if yes, continue (frequency = 3)
				shift_right(&temp, pattern_length);					// shift the copy of sequence by pattern_length bits
				if (masked_eq(&temp.m_lo, &seq512->m_lo, &mask))    // check if they match; if yes, return (frequency = 4)
					return 4;
				if (3 > curl) {
					if (last_curl_is_2)
						return 3;
					curl = 3;
				}
			}
			else if (2 > curl) {
				curl = 2;
			}
		}
	}
	return curl;													// no repetition found
}

// construct the sequence ( < 128 bits generator)
// __m256i* seq256 - generator of the sequence
//   uint8_t* len2 - length of the sequence
bool build_sequence128(__m256i* seq256, uint16_t* len2, uint8_t last_curl = 3) {
	//uint8_t last_curl = 3;										// if the curl would be > 3, we can break, thus we set the last curl to 3 in advance
	__m128i seq128 = _mm256_castsi256_si128(*seq256);				// transfer uint64_t to 128-bit register for larger lengths
	static const __m128i mask1_128 = _mm_set_epi64x(0, 1);
	for (; *len2 < 128; ++ * len2) {								// while within range of the register (256 bits)
		uint8_t curl = krul128(*len2, last_curl == 2, &seq128);		// find curl of the current sequence with cache
		if (curl == 1) {
			*seq256 = _mm256_set_m128i(_mm_setzero_si128(), seq128);
			return false;											// return the current length of the sequence if we encounter 1
		}
		if (curl == 4) {
			*seq256 = _mm256_set_m128i(_mm_setzero_si128(), seq128);
			return true;											// return the current length of the sequence + 1 if we encounter curl > 3, because the next curl will be 1
		}
		//shift_left_one(&seq128);									// shift the sequence left one bit
		shift_left(&seq128, 1);										// shift the sequence left one bit
		if (curl - 2)
			seq128 = _mm_or_si128(seq128, mask1_128);				// and add 1 if the curl is 3 (binary twin)

		last_curl = curl;											// save the curl
	}

	__m256i seq = _mm256_set_m128i(_mm_setzero_si128(), seq128);	// transfer 128-bit register to 256-bit register for larger lengths
	static const __m256i mask1 = _mm256_set_epi64x(0, 0, 0, 1);
	for (; *len2 < 255; ++ * len2) {								// while within range of the register (256 bits)
		uint8_t curl = krul256(*len2, last_curl == 2, &seq);		// find curl of the current sequence with cache
		if (curl == 1) {
			*seq256 = seq;
			return false;											// return the current length of the sequence if we encounter 1
		}
		if (curl == 4) {
			*seq256 = seq;
			return true;											// return the current length of the sequence + 1 if we encounter curl > 3, because the next curl will be 1
		}
		shift_left_one(&seq);										// shift the sequence left one bit
		if (curl - 2)
			seq = _mm256_or_si256(seq, mask1);						// and add 1 if the curl is 3 (binary twin)

		last_curl = curl;											// save the curl
	}

	// go to 512 bits
	m512 seq512(seq);												// construct 215-bit struct
	for (; *len2 < 512; ++ * len2) {								// while within range of the register (256 bits)
		uint8_t curl = krul512(*len2, last_curl == 2, &seq512);		// find curl of the current sequence with cache
		if (curl == 1) {
			*seq256 = seq512.m_lo; // TODO
			return false;											// return the current length of the sequence if we encounter 1
		}
		if (curl == 4) {
			*seq256 = seq512.m_lo; // TODO
			return true;											// return the current length of the sequence + 1 if we encounter curl > 3, because the next curl will be 1
		}
		push(&seq512, curl == 3);									// shift the sequence left one bit

		last_curl = curl;											// save the curl
	}

	std::cout << "encountered length 512" << std::endl;
	return false;
}

// construct the sequence ( < 64 bits generator)
// __m256i* seq256 - generator of the sequence
//   uint8_t* len2 - length of the sequence
bool build_sequence64(__m256i* seq256, uint16_t* len2) {
	uint8_t last_curl = 3;											// if the curl would be > 3, we can break, thus we set the last curl to 3 in advance
	uint64_t seq64 = _mm256_extract_epi64(*seq256, 0);
	// short curl for <= cache_bits
	for (; *len2 <= cache_bits; ++ * len2) {						// while within range of the 64-bit register
		uint8_t curl = cache[*len2][seq64 & mask64[*len2]];			// find curl of the current sequence with cache
		if (curl == 1) {
			*seq256 = _mm256_set_epi64x(0, 0, 0, seq64);
			return false;											// return the current length of the sequence if we encounter 1
		}
		if (curl == 4) {
			*seq256 = _mm256_set_epi64x(0, 0, 0, seq64);
			return true;											// return the current length of the sequence + 1 if we encounter curl > 3, because the next curl will be 1
		}
		seq64 <<= 1;
		if (curl - 2)
			seq64++;												// and add 1 if the curl is 3 (binary twin)

		last_curl = curl;											// save the curl
	}
	// short shift for <64
	for (; *len2 < 64; ++ * len2) {									// while within range of the 64-bit register
		uint8_t curl = krul64(*len2, last_curl == 2, seq64);		// find curl of the current sequence with cache
		if (curl == 1) {
			*seq256 = _mm256_set_epi64x(0, 0, 0, seq64);
			return false;											// return the current length of the sequence if we encounter 1
		}
		if (curl == 4) {
			*seq256 = _mm256_set_epi64x(0, 0, 0, seq64);
			return true;											// return the current length of the sequence + 1 if we encounter curl > 3, because the next curl will be 1
		}
		seq64 <<= 1;
		if (curl - 2)
			seq64++;												// and add 1 if the curl is 3 (binary twin)

		last_curl = curl;											// save the curl
	}

	*seq256 = _mm256_set_epi64x(0, 0, 0, seq64);
	return build_sequence128(seq256, len2, last_curl);
}

#pragma endregion

#pragma region MIRROR ALGORITHM

// build the sequence, calculate the tail length and potentially update best known value for tail length
void buildtail(__m256i* seq256, const uint8_t len, uint16_t* tail_length, bool* last_curl_is_4) {
	uint16_t len2 = len;											// len2 will become the length of the final sequence
	if (len <= 64)
		*last_curl_is_4 = build_sequence64(seq256, &len2);		// build this sequence (and check whether it terminates due to a 1 or a 4)
	else
		*last_curl_is_4 = build_sequence128(seq256, &len2);
	*tail_length = len2 - len;									// calculate tail length (without the final 1 and 4)
	if (*tail_length + *last_curl_is_4 > tails[len])			// if tail length (including a may-be-4) is higher than best known value of this thread, 
		tails[len] = *tail_length + *last_curl_is_4;			// update it
}

// returns whether the generator repeats itself (and could thus potentially profit from adding a prefix)
//__forceinline
bool checkperiod(__m256i* seq256, const uint8_t len, const uint8_t period) {
	__m256i temp = *seq256;										// copy of original sequence
	shift_right(&temp, period);									// shift the copy of sequence by period bits
	__m256i diff = _mm256_xor_si256(temp, *seq256);				// check for differing bits
	return (_mm256_testz_si256(diff, mask256[len - period]));					// and return whether there were any differing bits
}

// try to find valid extensions for this generator, and check their tails
void extend(const __m128i seq128, const uint8_t len, std::set<OptLen, compare_OptLen>& todo, const std::set<OptLen, compare_OptLen>& done) {
	__m256i seq256 = _mm256_set_m128i(_mm_setzero_si128(), seq128);			// set seq64 in 256-bit register
	bool last_curl_is_4;										// will hold whether the sequence terminated due to a 1 or a 4
	uint16_t tail_length = 0;									// will hold the tail length of the sequence
	buildtail(&seq256, len, &tail_length, &last_curl_is_4);		// build tail, calculate tail length, set last_curl_is_4 (and update 'tails' if new record found)
	if (len < length) {											// if we didn't reach max length, let's try to add a prefix to the generator
		for (uint8_t i = 0; i < tail_length + !last_curl_is_4; ++i) {
			int pos = tail_length - i;												// position of the bit of interest
			uint8_t k = 3;															// if the curl at this position was 2, we can improve to it to 3
			if ((i == (tail_length - 1 + !last_curl_is_4)) && (!last_curl_is_4)) {	// if this position is the last element and it was a 1: k = 2
				k = 2;
			}
			else if (seq256.m256i_u64[(pos - 1) / 64] & (1ull << ((pos - 1) % 64)))	// and if the curl at this position was a 3, skip,
				continue;															// because improving to curl = 4 is useless
			__m256i copy = seq256;							// copy of sequence up to pos bits
			shift_right(&copy, pos);							// copy of sequence up to pos bits
			uint8_t lower_bound = ((i + len) / k) + 1;								// smallest possible period
			uint8_t upper_bound = ((i + length) / k);								// largest possible period
			for (uint8_t p = lower_bound; p <= upper_bound; ++p) {
				if (checkperiod(&copy, i + len, p)) {								// check if the copy repeats itself
					uint8_t modulo = (i + len) % p;									// modulo is that part of the repetition that is already present
					__m256i temp = copy;					// so we can discard these bits
					shift_right(&temp, modulo);					// so we can discard these bits
					__m128i prefix = _mm256_castsi256_si128(temp);
					prefix = _mm_and_si128(prefix, mask128[p - modulo]);			// and we select that part of the repetition that needs to be prefixed
					shift_left(&prefix, len);									// TODO: SHIFT > 64 BITS
					__m128i option = _mm_or_si128(prefix, seq128);					// and after prefixing, this is the next option
					uint8_t opt_len = len + p - (modulo);							// and we calculate its length
					OptLen s(option, opt_len);
					auto done_it = done.find(s);									// now, we try to find this option in the sequences we already did
					if (done_it == done.end()) {									// and if "opt" is NOT already processed
						todo.insert(s);												// we add it to the todo (and if it was already there, we discard it)
					}
				}
			}
		}
	}
}

void find_possible_records(int id) {
	std::set<OptLen, compare_OptLen> todo, done;
	uint64_t start = (1ull << (length / fraction)) * id / thread_count;		// start value for this thread (assuming we start from id = 0)
	uint64_t stop = (1ull << (length / fraction)) * (id + 1) / thread_count;	//  stop value for this thread (assuming we start from id = 0)
	for (uint64_t seq64 = start; seq64 < stop; ++seq64) {
		__m128i seq128 = _mm_set_epi64x(0, seq64);
		todo.insert(OptLen(seq128, length / fraction));
		while (!todo.empty()) {
			OptLen s = *todo.begin();					// get the first option/length struct
			done.insert(s);								// put that struct into "done"
			todo.erase(todo.begin());					// remove that from "todo"
			extend(s.opt, s.len, todo, done);
		}
		done.clear();
	}
	{	// at the end, we update the global array of tails (under a thread lock) if tails in these threads were higher
		const std::lock_guard<std::mutex> l(m_tails);
		for (int i = 0; i <= length; ++i)
			global_tails[i] = std::max(global_tails[i], tails[i]);
	}
}

// starts each thread with a unique ID and joins them
void multi_threader() {
	std::vector<std::thread> thread_vector;
	for (int id = 0; id < thread_count; ++id)
		thread_vector.emplace_back(std::thread(find_possible_records, id));
	for (auto& th : thread_vector) th.join();
}
#pragma endregion

#pragma region CACHE BUILDING
void update_curls(int m, std::array<uint8_t, 8>& curls, uint8_t curl) {
	unsigned long index = 0, base = 0;
	while (_BitScanForward(&index, m)) {
		if (curls[index + base] < curl)
			curls[index + base] = curl; // set corresponding curls, if they are less
		base += ++index;
		m >>= index;
	}
}

void update_curls4(int m, std::array<uint8_t, 8>& curls) {
	unsigned long index = 0, base = 0;
	while (_BitScanForward(&index, m)) {
		curls[index + base] = 4; // set corresponding curls, if they are less
		base += ++index;
		m >>= index;
	}
}

// build 32-bit cache
void krul_for_cache(uint32_t seq16, uint32_t cache_width) {
	static const __m256i mask_one = _mm256_set_epi32(1, 1, 1, 1, 1, 1, 1, 1);
	__m256i mask = mask_one;
	std::array<uint8_t, 8> curls = { 1, 1, 1, 1, 1, 1, 1, 1 };
	__m256i seq256 = _mm256_set_epi32(seq16 + 7, seq16 + 6, seq16 + 5, seq16 + 4, seq16 + 3, seq16 + 2, seq16 + 1, seq16);
	int len = 1;
	int limit2 = cache_width / 2;
	int limit3 = cache_width / 3;
	int limit4 = cache_width / 4;
	int m2 = 0, m3 = 0, m4 = 0;
	int m2done = 0, m3done = 0, m4done = 0;
	// do <= 1/2 of cache_width, for possible freq 
	for (; len <= limit2; ++len) {
		__m256i target = _mm256_and_si256(seq256, mask);
		__m256i temp = _mm256_srli_epi32(seq256, len);
		__m256i diff = _mm256_and_si256(temp, mask);
		diff = _mm256_cmpeq_epi32(diff, target);
		m2 = _mm256_movemask_ps(_mm256_castsi256_ps(diff));
		if (m2) { // any of the tails are equal? keep going...
			if (len <= limit3) { // can we do 3?
				temp = _mm256_srli_epi32(temp, len);
				diff = _mm256_and_si256(temp, mask);
				diff = _mm256_cmpeq_epi32(diff, target);
				m3 = _mm256_movemask_ps(_mm256_castsi256_ps(diff));
				m3 &= m2; // only consider if already found 2
				if (m3) { // any of the tails are equal? keep going...
					if (len <= limit4) { // can we do 4?
						temp = _mm256_srli_epi32(temp, len);
						diff = _mm256_and_si256(temp, mask);
						diff = _mm256_cmpeq_epi32(diff, target);
						m4 = _mm256_movemask_ps(_mm256_castsi256_ps(diff));
						m4 &= m3; // only consider if already found 3
						m4 &= ~m4done; // exclude already marked 4
						if (m4) { // any of the tails are equal? found 4
							update_curls4(m4, curls);
							m4done |= m4;
						}
						m3 &= ~m4; // exclude m4 from m3
					}
					m3 &= ~m3done; // exclude already marked 3
					update_curls(m3, curls, 3);
					m3done |= m3;
				}
				m2 &= ~m3; // exclude m3 from m2
			}
			m2 &= ~m2done; // exclude already marked 2
			update_curls(m2, curls, 2);
			m2done |= m2;
		}
		mask = _mm256_slli_epi32(mask, 1); // shift left by 1
		mask = _mm256_or_si256(mask, mask_one); // add 1
	}

	uint8_t* pcurl = &cache[cache_width][seq16];
	for (int i = 0; i < 8; ++i) {
		*pcurl++ = curls[i];
	}
}

void build_cache_thread(uint32_t start, uint32_t stop, uint32_t cache_width)
{
	for (uint32_t seq64 = start; seq64 < stop; seq64 += 8) // do 8 32-bit at a time
		krul_for_cache(seq64, cache_width);
}

void build_cache() {
	// single thread for very short cache because of the overhead
	for (int cache_width = 1; cache_width <= 16; ++cache_width) {
		cache[cache_width] = (uint8_t*)malloc(std::max(8, (1 << cache_width))); // at least 8 bytes
		const uint32_t end = 1 << cache_width;
		for (uint32_t seq16 = 0; seq16 < end; seq16 += 8) {
			krul_for_cache(seq16, cache_width);
		}
	}
	// multi-threading
	for (uint32_t cache_width = 17; cache_width <= cache_bits; ++cache_width) {
		cache[cache_width] = (uint8_t*)malloc(1ull << (cache_width));

		const uint32_t begin = 0;
		const uint32_t end = 1ull << cache_width;
		const uint32_t span = ((end - begin) / thread_count) & ~0x7; // make sure it's divisible by 8

		uint32_t start = begin;	      // start value for this thread
		uint32_t stop = start + span; // end value for this thread

		std::vector<std::thread> thread_vector;

		// start the threads and join them
		for (int i = 0; i < thread_count; ++i) {
			if (i == thread_count - 1)
				stop = end; // avoid rounding down error
			thread_vector.emplace_back(std::thread(build_cache_thread, start, stop, cache_width));
			start += span;
			stop += span;
		}
		for (auto& th : thread_vector) { th.join(); }
	}
	cache30 = cache[cache_bits];
}
#pragma endregion

int main() {
	//std::ofstream file("timing.txt", std::ios::app);
	//auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	////file << "\nStart: " << ctime(&timenow);
	//file << "Length: " << length << ", cache: " << cache_bits << ", threads: " << thread_count << std::endl;
	//file.close();

	const std::vector<int> tails_article = { 0,  // extra 0 for easy access
		  0,   2,   2,   4,   4,   8,   8,  58,  59,  60, 112, 112, 112, 118, 118, 118, 118, 118, 119, 119, // 20
		119, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, // 40
		120, 120, 120, 120, 120, 120, 120, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, // 60
		131, 131, 131, 131, 131, 131, 131, 132,	132, 132, 132, 132, 133, 133, 133, 133, 173, 173, 173, 173, // 80	73-75: Steven
		173, 173, 173, 173, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, // 100	Steven
		179, 179, 179, 179, 179, 179, 179,																	// 108	Steven
	}; 
	auto t1 = std::chrono::high_resolution_clock::now();
	build_cache();
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Build cache: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msec" << std::endl;
	multi_threader();
	//find_possible_records(0);
	auto t3 = std::chrono::high_resolution_clock::now();
	for (int i = length / fraction; i <= length; i++) {
		if (global_tails[i] < global_tails[i - 1])
			global_tails[i] = global_tails[i - 1];
		if (i < (int)tails_article.size()) {
			if (global_tails[i] != tails_article[i])
				std::cout << "WRONG: " << i << ", " << global_tails[i] << std::endl;
		}
		else
			std::cout << "NEW: " << i << ", " << global_tails[i] << std::endl;
	}
	std::cout << "[" << length << "] : " <<
		std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << " msec" << std::endl;

	//std::ofstream file("timing.txt", std::ios::app);
	//auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	////file << "Finish: " << ctime(&timenow);
	//file << "Length: " << length << ", cache: " << cache_bits << ", tail: " << tail << ", time: " <<
	//    std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count() << " sec" << std::endl;
}
