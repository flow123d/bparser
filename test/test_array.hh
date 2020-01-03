/*
 * test_array.hh
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */


#include <chrono>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <random>
#include "vec.hh"

/**
 * Question one
 */


/**
Timing (minimum over 3 runs):

 ** -O3 **
Sum, subarray, masked, index access.
  Time: 13 Sum: 2.49231e+09
Sum, subarray, indirect, index access.
  Time: 6 Sum: 2.49231e+09
Sum, subarray, indirect, ptr loop.
  Time: 6 Sum: 2.49231e+09
Sum, subarray, indirect, iterator loop.
  Time: 6 Sum: 2.49231e+09
Sum, full array, direct, foreach.
  Time: 6 Sum: 4.995e+09
Sum, const array, indirect, index access.
  Time: 5 Sum: 0
Sum, Array, sumloop, subarray.
  Time: 6 Sum: 2.49231e+09
Sum, Array, sumloop, full.
  Time: 6 Sum: 124750
Sum, Array, sumloop, const.
  Time: 5 Sum: 0


** -O2 **
Sum, subarray, masked, index access.
  Time: 15 Sum: 2.49231e+09
Sum, subarray, indirect, index access.
  Time: 14 Sum: 2.49231e+09
Sum, subarray, indirect, ptr loop.
  Time: 14 Sum: 2.49231e+09
Sum, subarray, indirect, iterator loop.
  Time: 15 Sum: 2.49231e+09
Sum, full array, direct, foreach.
  Time: 16.5 Sum: 4.995e+09
Sum, const array, indirect, index access.
  Time: 15 Sum: 0
Sum, Array, sumloop, subarray.
  Time: 6 Sum: 2.49231e+09
Sum, Array, sumloop, full.
  Time: 13 Sum: 124750
Sum, Array, sumloop, const.
  Time: 13 Sum: 0

Conclusion: When using -O3 all subset loop implementations takes about a same time. For -O2 the indirect access is probably
more costly. Using a subset pointer seems to be good choice to optimize memory consumption of continuous and constant arrays
while keeping the implementation assume for all three variants and doesn't be suboptimal e.g. for the constant case;
The resulting times corresponds to about 10e10 flops per 6 seconds:
performace: 1.6 GFlops
 on single core running at 2.7GHz,
theoretical peak performance:  2.7 GFlops

Further notes:
- Resulting speed is extremally sensitive, using SSE and just little bit different code results in the
"-O2" behavior, i.e. just the case "Sum, Array, sumloop, subarray." gives the short time.
Further investigation of the assembly is necessary. Possible detailed optimization can leed to even better times.


SSE, AVX notes:
- are not applicable with indirect adressing, the data must be in continuous blocks (2 or 4 doubles)
- SSE can compute with 2 doubles or 4 floats
- AVX2 can compute with 4 doubles or 8 floats
- AVX512 even with 8 doubles ...

Conclusion, if we succeed to store data into continuous blocks of 2,4, or 8 doubles,
 we can achieve significant speedup. See this:
> c++ -pedantic-errors -Wall -Wextra -Werror -Wno-long-long -O3 -std=c++11 -mavx2 -ffast-math -ftree-vectorizer-verbose=10 -I include test_array_loop.cc
> ./a.out
Sum, subarray, masked, index access.
  Time: 15 Sum: 2.49231e+09
Sum, subarray, indirect, index access.
  Time: 14 Sum: 2.49231e+09
Sum, subarray, indirect, ptr loop.
  Time: 14 Sum: 2.49231e+09
Sum, subarray, indirect, iterator loop.
  Time: 15 Sum: 2.49231e+09
!!! Sum, full array, direct, foreach.
  Time: 1.5 Sum: 4.995e+09
!!! Sum, const array, indirect, index access.
  Time: 1 Sum: 0s
Sum, Array, sumloop, subarray.
  Time: 6 Sum: 2.49231e+09
Sum, Array, sumloop, full.
  Time: 6 Sum: 124750
Sum, Array, sumloop, const.
  Time: 6 Sum: 0

For flat and cost explicit implementation (i.e. without indirect addrassing)
we get about: 10 GFlops !!

We can get close to it with Array having continuous simd blocks:

Sum, Array, sumloop, subarray.
  Time: 1 Sum: 2.61776e+09
Sum, Array, sumloop, full.
  Time: 3 Sum: 130816
Sum, Array, sumloop, const.
  Time: 3 Sum: 768


Conclusion: We need at least small continuous blocks. Not clear why full and constant is slower.

TODO:
- implement Array with all functions necessary in exprtk
- for individual operations or small expressions make comparison:
  - against hand optimized code
  - count number of flops and compare against theoretical power
  - compare both on small vectors (that fit into cache) and large vectors

**/

typedef unsigned int uint;
const uint simd_block_size=4;
const uint array_max_size = 1024;
const uint array_sub_size = 768;


#define FOR(i, j) \
	for(uint i=0; i < size_; ++i)  \
		for(uint j=0; j<simd_block_size; j++)

struct Array {
	Array(double *ptr, uint *subset, uint size)
	: ptr_(ptr),
	  subset_(subset),
	  size_(size)
	{}

	inline void reinit(double *ptr) {
		ptr_ = ptr;
	}

	inline double & at(uint i, uint j) {
		return *(ptr_ + *(subset_ + i) + j);
	}

	inline double sum_loop() {
		// one add, short vector: 1 3 3
		// one add, long vector: 6 3 4
		// one mult one add, long vector: 9 4 2
		double sum = 0;
		FOR(i,j) {
			sum += 4.1 * at(i,j) + 1.2 ;
		}
		return sum;
	}

	inline void add_loop(Array &other) {
		// one add, long vector: 9 9 9
		FOR(i,j) {
			other.at(i,j) = at(i,j) + at(i,j);
		}
	}

	double *ptr_;
	uint *subset_;
	uint size_;


};



template <class Expr>
struct ConstArray {
	ConstArray() {
		for(uint i=0; i < simd_block_size; ++i)
			a[i] = 3.14;
	}

	inline double run() {
		double sum = 0;
		for(uint i=0; i < simd_block_size; ++i)
			Expr
			sum += a[i];
		}

	}

	double a[simd_block_size];
};


template <class Expr>
struct FlatArray {
	FlatArray() {
		size = array_sub_size;
	}
	inline void run() {
		for(uint i=0; i < array_sub_size/simd_block_size; ++i)
			for(uint j=0; j < simd_block_size; ++j)
				Expr::expr(a, b, c[j], d[j]);
	}

	double check() {
		double sum = 0;
		for(uint i=0; i < size; ++i)
			sum += res[i]
	}

	uint size;
	double a1[array_max_size];
	double a2[array_max_size];
	double a3[array_max_size];
	double b[array_max_size];
	double c1[simd_block_size];
	double c2[simd_block_size];
	double c3[simd_block_size];
	double d[simd_block_size];
	double res[array_max_size];
};

struct SubsetArray {

};


struct MaskedArray {
	MaskedSum() {

	}

	inline std::string name() {
		return "Sum, subarray, masked, index access.";
	}

	inline uint n_flop() {
		return array_sub_size;
	}

	inline double run() {
		double sum = 0;
		for(uint i=0; i < size; ++i)
		  if (mask[i]) {
			sum += a[i];
		}

	}
};


struct SumExpr {
	inline void add(double &sum, double value) {
		sum += value;
	}
};

struct SumProdExpr {
	inline void add(double &sum, double value) {
		sum += value;
	}
};



// Make a real case subset: few continuous blocks of random size.
std::vector<Block> make_subset() {
	uint total_size = array_max_size;
	uint subset_size = array_sub_size;
	uint n_blocks = 7;

	std::default_random_engine generator(1011);
	double sum_sizes = 0;
	double sum_shifts = 0;
	std::vector<uint> sizes;
	std::vector<uint> shifts;
	for(uint i; i< n_blocks; ++i) {
		uint val = generator();
		sum_sizes += val;
		sizes.push_back(val);
		uint val = generator();
		sum_shifts += val;
		shifts.push_back(val);
	}

	std::vector<Block> subset;
	uint pos = 0;
	double total_shift = array_max_size - array_sub_size;
	uint actual_subset_size = 0;
	for(uint i; i< n_blocks; ++i) {
		uint shift = uint(total_shift * (shifts[i] / sum_shifts));
		uint start = pos + shift;
		uint size = uint(subset_size * (sizes[i] / sum_sizes));
		actual_subset_size += size;
		Block b = {start, size};
		subset.push_back(b);
		pos += shift + size;
	}
	uint diff = subset_size - actual_subset_size
	ASSERT(diff > 0);
	subset.back().size += diff;
}



struct ConstAssign {
	ConstAssign() {
		c = 3.14;

	}

	double c;

	double res[array_max_size];
};

template <class BenchCase>
void measure() {
	const uint n_repeat = 10000;
	BenchCase bc;

	std::cout << bc.name() << "\n";
	auto start_time = std::chrono::high_resolution_clock::now();

	double check_sum = 0;
	for(uint j=0; j<n_repeat; ++j) {
		check_sum += bc.run();
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double time  = (end_time - start_time)/std::chrono::milliseconds(1);
	std::cout << "  perf. [GFlop/s]: " << n_repeat * bc.n_flop / time << "\n";
	std::cout << "  time [ms]:       " << time << "\n";
	std::cout << "  check sum:       " << check_sum << " << \n";

}



int main() {

	uint size = 1024;

	std::vector<bool> mask(size, false);
	for(uint i=0; i < size; ++i) mask[i] = bool(std::rand() % 2);

	std::vector<uint> subset;
	std::vector<uint> full_subset;
	std::vector<uint> const_subset;

	for(uint i=0; i < size; ++i) {
		if (mask[i]) subset.push_back(i);
		full_subset.push_back(i);
		const_subset.push_back(0);
	}

	std::vector<double> a(n_repeat * size);
	for(uint i=0; i < size; ++i) a[i] = i % size;

	{
		std::cout << "Sum, subarray, masked, index access." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			for(uint i=0; i < size; ++i)
			  if (mask[i]) {
				sum += a[i];
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, subarray, indirect, index access." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			for(uint i : subset) {
				sum += a[i];
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, subarray, indirect, ptr loop." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		double *a_ptr = &(a[0]);
		for(uint j=0; j<n_repeat; ++j) {
			uint *subset_ptr = &(*(subset.begin()));
			uint *subset_end = subset_ptr + subset.size();
			for(; subset_ptr != subset_end; ++subset_ptr) {
				sum += *(a_ptr + *subset_ptr);
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, subarray, indirect, iterator loop." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			double *a_ptr = &(a[0]);
			for(auto it=subset.begin(); it != subset.end(); ++it) {
				sum += *(a_ptr + *it);
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, full array, direct, foreach." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			double *a_ptr = &(a[0]);
			for(uint i=0; i<subset.size(); ++i) {
				sum += *(a_ptr + i);
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time/2 << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, const array, indirect, index access." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			double *a_ptr = &(a[0]);
			for(uint i=0; i<subset.size(); ++i) {
				sum += *a_ptr;
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	// Make simd blocks
	subset.clear();
	full_subset.clear();
	const_subset.clear();

	for(uint i=0; i < size; i+=simd_block_size) {
		if (mask[i]) subset.push_back(i);
		full_subset.push_back(i);
		const_subset.push_back(0);
	}


	std::vector<Array> workspace = {
			Array(&(a[0]), &(*(subset.begin())), subset.size()),
			Array(&(a[0]), &(*(full_subset.begin())), full_subset.size()/2),
			Array(&(a[0]), &(*(const_subset.begin())), const_subset.size()/2)};


	{
		std::cout << "Sum, Array, sumloop, subarray." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			workspace[0].reinit(&(a[j*size]));
			sum += workspace[0].sum_loop();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, Array, sumloop, full." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			workspace[1].reinit(&(a[j*size]));
			sum += workspace[1].sum_loop();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}

	{
		std::cout << "Sum, Array, sumloop, const." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		double sum = 0;
		for(uint j=0; j<n_repeat; ++j) {
			workspace[2].reinit(&(a[j*size]));
			sum += workspace[2].sum_loop();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << " Sum: " << sum << "\n";
	}


	std::vector<double> b(n_repeat*size);
	Array b_array(&(b[0]), &(*(subset.begin())), subset.size());

	{
		std::cout << "Sum, Array, addloop, subarray." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint j=0; j<n_repeat; ++j) {
			b_array.reinit(&(b[j*size]));
			workspace[0].add_loop(b_array);
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time <<  "\n";
	}

	{
		std::cout << "Sum, Array, addloop iter, full." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint j=0; j<n_repeat; ++j) {
			b_array.reinit(&(b[j*size]));
			workspace[1].add_loop(b_array);
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time <<  "\n";
	}

	{
		std::cout << "Sum, Array, addloop iter, const." << "\n";
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint j=0; j<n_repeat; ++j) {
			b_array.reinit(&(b[j*size]));
			workspace[2].add_loop(b_array);
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		double time  = (end_time - start_time)/std::chrono::milliseconds(1);
		std::cout << "  Time: " << time << "\n";
	}

}


