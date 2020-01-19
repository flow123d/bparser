/*
 * test_design.cc
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */

/**
 * Observation:
 * - calculation on flat array is faster then on subset array,
 *   compiler have problem to use vectorization in the subset case
 * - However due to copies for the flat variant, this performance drawback is pa
 *
 * - direct calculation is much fater, fastest for complex expressions, where there is just single load of values into
 *   AVX unit and then multiple operations is performed.
 *
 * - for complex expressions even three times faster
 * - problem is that all loops are condensed into a single one, need to compute on different data
 * - examination of the assembly shows that the subset variant is unable to vectorize, try  a loop modification.
 * - then it is curious that the flat variant is only two times faster not 4 times.
 *
 * TODO:
 */

#include <chrono>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <array>
#include <malloc.h>
#include "assert.hh"

typedef unsigned int uint;
const uint simd_block_size=4;

const uint array_max_size = 1024 / simd_block_size;
const uint array_sub_size = 768 / simd_block_size;

//const uint array_max_size = 16 / simd_block_size;
//const uint array_sub_size = 8 / simd_block_size;

typedef double double4 __attribute__((__vector_size__(32)));
//typedef std::array<double, simd_block_size> double4;





template <class BenchCase>
void measure(BenchCase &bc) {

	const uint n_repeat = 1000000;
	//const uint n_repeat = 1;


	auto start_time = std::chrono::high_resolution_clock::now();
	for(uint j=0; j<n_repeat; ++j) {
		uint shift = j % bc.d.n_arrays;
		bc.run_array(shift);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	double time_array  = (end_time - start_time)/std::chrono::milliseconds(1);

	start_time = std::chrono::high_resolution_clock::now();
	for(uint j=0; j<n_repeat; ++j) {
		uint shift = j % bc.d.n_arrays;
		bc.run_direct(shift);
	}
	end_time = std::chrono::high_resolution_clock::now();
	double time_direct  = (end_time - start_time)/std::chrono::milliseconds(1);

	double check_sum = bc.d.check();

	// Full output
//	std::cout << typeid(BenchCase).name() << "\n";
//	std::cout << "  perf. [GFlop/s]: " << n_repeat * bc.n_flop * array_sub_size / time_array / 1000000 << "\n";
//	std::cout << "  (direct perf)  : " << n_repeat * bc.n_flop * array_sub_size / time_direct / 1000000 << "\n";
//	std::cout << "  time [ms]      : " << time_array << "\n";
//	std::cout << "  time direct    : " << time_direct << "\n";
//	std::cout << "  check sum:       " << check_sum << "\n";

	// condensed
	std::cout << typeid(BenchCase).name() << "  perf. [GFlop/s]: " << n_repeat * bc.n_flop * array_sub_size * simd_block_size / time_array / 1000000 << "\n";
	std::cout << typeid(BenchCase).name() << "  (direct perf.) : " << n_repeat * bc.n_flop * array_sub_size * simd_block_size / time_direct / 1000000 << "\n";

}


/**
 * Question one: Compare speed of calculation on subset vector and calculation on flat vector.
 */
void make_subset(uint seed, uint *subset) {
	srand(seed);
	uint n_blocks = array_sub_size;
	uint shift = 0;
	uint total_shift = (array_max_size - array_sub_size);
	for(uint i_block = 0; i_block < n_blocks; ++i_block) {
		int add_shift = int(total_shift * (rand() % array_max_size) / array_max_size);
		shift += add_shift;
		total_shift -= add_shift;
		subset[i_block] = shift; // alligned shifts
		shift += 1;
		//std::cout << i_block << " " << subset[i_block] << " " << total_shift << "\n";
		BP_ASSERT(subset[i_block] < array_max_size );
	}
}

#define SubsetFor \
		for(uint i=0; i<n_blocks;  ++i) \
			for(uint j=0; j<simd_block_size; ++j)


struct SubsetArray {
	typedef SubsetArray Array;

	SubsetArray() {
	}


	void input(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
		subset = subset_ptr;
		values = ptr;
	}

	void output(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
		subset = subset_ptr;
		values = ptr;
	}

	void middle(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
		subset = subset_ptr;
		values = ptr;
	}

	void copy_out(double4 *ptr, uint *subset_ptr, uint subset_size) {
//		for(uint i=0; i<n_blocks; ++i)
//			for(uint j = 0, idx=0; idx = subset[i] + j, j < simd_block_size; ++j) {
//				std::cout << i << " " << j << " " << values[idx] << "\n";
//			}

	}

	uint idx(uint i) const
	{ return subset[i];}

	void print() {
		SubsetFor {
				std::cout << idx(i) << " " << j << " " << values[idx(i)][j] << "\n";
			}
	}

	double sum() {
		sum_[0] = sum_[1] = sum_[2] = sum_[3] = 0.0;
		SubsetFor {
			sum_[j] += values[idx(i)][j];
		}

		return sum_[0] + sum_[1] + sum_[2] + sum_[3];
	}

	void set_sum(const Array &a, const Array &b) {
		SubsetFor {
			values[idx(i)][j] = a.values[a.idx(i)][j] + b.values[b.idx(i)][j];
		}

		// not vectorized
	}

	void add(const Array &a) {
		SubsetFor
		    values[idx(i)][j] += a.values[a.idx(i)][j];
	}

	void scalar_add(double4 cc) {
		SubsetFor
			values[idx(i)][j] += cc[j];
		// not vectorized
	}

	void mult(const Array &a, const Array &b) {
		SubsetFor
			values[idx(i)][j] = a.values[a.idx(i)][j] * b.values[b.idx(i)][j];
	}


	void scale(const Array &a) {
		SubsetFor
			values[idx(i)][j] *= a.values[a.idx(i)][j];
	}

	void scalar_scale(double4 cc) {
		SubsetFor
			values[idx(i)][j] *= cc[j];
	}

	void sqrt_(const Array &a) {
		SubsetFor
			values[idx(i)][j] = sqrt(a.values[a.idx(i)][j]);
	}

	double4 sum_;
	uint n_blocks;
	uint *subset;
	double4 * values;
};

#define FlatFor \
		for(uint idx=0; idx<n_blocks; idx+=1) \
			for(uint j=0; j<simd_block_size; ++j)


struct FlatArray {
	typedef FlatArray Array;

	FlatArray () {
		// padding the sub size
		n_blocks = (array_sub_size + simd_block_size -1) / simd_block_size;
		//Make own storage
	}

	void input(double4 *ptr, uint *subset_ptr, uint subset_size) {
		// Copy
		n_blocks = subset_size;
		for(uint i=0; i<n_blocks; ++i) {
			uint f_idx = i;
			uint s_idx = subset_ptr[i];
			for(uint j=0; j<simd_block_size; ++j) {
				values[f_idx][j] = ptr[s_idx][j];
			}
				//std::cout << "copy: " << f_idx << " <- " << s_idx <<  " " << ptr[s_idx] << "\n";
		}

	}

	void output(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
	}

	void middle(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
	}


	void copy_out(double4 *ptr, uint *subset_ptr, uint subset_size) {
		n_blocks = subset_size;
		for(uint i=0; i<n_blocks; ++i) {
			uint f_idx = i;
			uint s_idx = subset_ptr[i];
			for(uint j=0; j<simd_block_size; ++j) {
				ptr[s_idx][j] = values[f_idx][j];
			}
			// std::cout << i << " " << j << " " << values[f_idx] << "\n";

		}
	}

	void print() {
		n_blocks = array_sub_size / simd_block_size;
		for(uint i=0; i<n_blocks; ++i) {
			for(uint j = 0, idx=0; j < simd_block_size; ++j) {
				std::cout << i << " " << j << " " << values[i][j] << "\n";
			}
		}
	}

	double sum() {
		sum_[0] = sum_[1] = sum_[2] = sum_[3] = 0.0;
		FlatFor {
			//std::cout << i << " " << idx << " " << values[idx] << "\n";
			sum_[j] += values[idx][j];
		}

		return sum_[0] + sum_[1] + sum_[2] + sum_[3];
	}

	void set_sum(const Array &a, const Array &b) {
		FlatFor
		    values[idx][j] = a.values[idx][j] + b.values[idx][j];
		// assembly: vectorized
	}

	void add(const Array &a) {
		FlatFor
			values[idx][j] += a.values[idx][j];
	}

	void scalar_add(double4 a) {
		FlatFor
			values[idx][j] += a[j];
	}

	void mult(const Array &a, const Array &b) {
		FlatFor
			values[idx][j] = a.values[idx][j] * b.values[idx][j];
	}


	void scale(const Array &a) {
		FlatFor
			values[idx][j] *= a.values[idx][j];
	}

	void scalar_scale(double4 b) {
		FlatFor
			values[idx][j] *= b[j];
	}

	void sqrt_(const Array &a) {
		FlatFor
			values[idx][j] = sqrt(a.values[idx][j]);
	}


	double4 values[array_max_size];
	double4 sum_;
	uint n_blocks;
};


template <class A>
struct Expr {
	static const uint n_arrays = 5;
	static const uint n_temp = 4;
	static const uint n_blocks = array_sub_size;

	Expr() {
		//std::cout << "Expr constr" << "\n";
		base_size = (2 * n_arrays + n_temp) * array_max_size;
		base = (double4 *)memalign(sizeof(double) * simd_block_size, sizeof(double) * simd_block_size * base_size);
		BP_ASSERT(base !=nullptr);
		// Initialization
		for(uint i=0; i<base_size; ++i)
			for(uint j = 0, idx=0; j < simd_block_size; ++j) {
			base[i][j] = i * simd_block_size + j;
		}

		for(uint i = 0; i < 2 * n_arrays; ++i) {
			a_ptr[i] = base + i * array_max_size;
		}

		for(uint j=0;j<n_temp;++j) {
			double cc = 3.14 * j;
			for(uint k=0; k<simd_block_size; ++k)
				constants[j][k] = cc;
		}

		for(uint j=0;j<4;++j)
			temp[j] = base + (2 * n_arrays + j) * array_max_size;

		uint seed = 1234;

		make_subset(seed, subset);
		for(uint i=0; i<n_blocks; ++i)
			flat_subset[i] = i;
	}

	double check() {
		sum_[0] = sum_[1] = sum_[2] = sum_[3] = 0.0;
		for(uint i=0; i < 2 * n_arrays  * array_max_size; ++i)
			for(uint j = 0, idx=0; j < simd_block_size; ++j) {
			//std::cout << i << " " << idx << " " << values[idx] << "\n";
			sum_[j] += base[i][j];
		}
		return sum_[0] + sum_[1] + sum_[2] + sum_[3];
	}

	~Expr() {
		free(base);
	}

	double4 sum_;
	double4 constants[n_temp];
	A a[n_arrays];
	double4 * temp[n_temp];
	double4 * a_ptr[2 * n_arrays];
	double4 * base;
	uint base_size;
	uint subset[n_blocks];
	uint flat_subset[n_blocks];
};



template <class A>
struct Sum {
	Expr<A> d;
	static const uint n_flop = 1;

	void run_array(uint shift) {

		d.a[0].output(d.a_ptr[shift+0], d.subset, d.n_blocks);
		d.a[1].input(d.a_ptr[shift+1], d.subset, d.n_blocks);
		d.a[2].input(d.a_ptr[shift+2], d.subset, d.n_blocks);

		d.a[0].set_sum( d.a[1], d.a[2]);

		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);
	}

	void run_direct(uint shift) {
		//std::cout << "direct" << std::endl;
		double *a[3];
		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
		a[1] = (double *)(& (d.a_ptr[shift+1][0]));
		a[2] = (double *)(& (d.a_ptr[shift+2][0]));
		for(uint i=0; i<array_sub_size * simd_block_size; ++i) {
			a[0][i] = a[1][i] + a[2][i];
		}
	}

};


template <class A>
struct Add {
	Expr<A> d;
	static const uint n_flop = 1;

	void run_array(uint shift) {

		d.a[0].input(d.a_ptr[shift+0], d.subset, d.n_blocks);
		d.a[1].input(d.a_ptr[shift+1], d.subset, d.n_blocks);

		d.a[0].add( d.a[1]);

		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);

	}

	void run_direct(uint shift) {
		double *a[2];
		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
		a[1] = (double *)(& (d.a_ptr[shift+1][0]));
		for(uint i=0; i<array_sub_size * simd_block_size; ++i)
			a[0][i] += a[1][i];
	}

};





template <class A>
struct AddScalar {
	Expr<A> d;
	static const uint n_flop = 1;

	void run_array(uint shift) {

		d.a[0].input(d.a_ptr[shift+0], d.subset, d.n_blocks);

		d.a[0].scalar_add( d.constants[0]);

		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);
	}

	void run_direct(uint shift) {
		double *a[1];
		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
		for(uint i=0; i<array_sub_size * simd_block_size; ++i)
			a[0][i] += 3.14;
	}

};



template <class A>
struct NormSqr {
	Expr<A> d;
	static const uint n_flop = 5;

	void run_array(uint shift) {
		d.a[0].output(d.a_ptr[shift+0], d.subset, d.n_blocks);
		d.a[1].input(d.a_ptr[shift+1], d.subset, d.n_blocks);
		d.a[2].input(d.a_ptr[shift+2], d.subset, d.n_blocks);
		d.a[3].input(d.a_ptr[shift+3], d.subset, d.n_blocks);
		d.a[4].middle(d.temp[0], d.flat_subset, d.n_blocks);


		d.a[0].mult( d.a[1], d.a[1]);

		d.a[4].mult( d.a[2], d.a[2]);
		//d.a[4].print();
		d.a[0].add( d.a[4]);

		d.a[4].mult( d.a[3], d.a[3]);
		d.a[0].add( d.a[4]);


		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);


	}

	void run_direct(uint shift) {
		double *a[4];
		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
		a[1] = (double *)(& (d.a_ptr[shift+1][0]));
		a[2] = (double *)(& (d.a_ptr[shift+2][0]));
		a[3] = (double *)(& (d.a_ptr[shift+3][0]));
		for(uint i=0; i<array_sub_size * simd_block_size; ++i)
			a[0][i] = a[1][i]*a[1][i] + a[2][i]*a[2][i] + a[3][i]*a[3][i];
	}

};


template <class A>
struct Sqrt {
	Expr<A> d;
	static const uint n_flop = 1;

	void run_array(uint shift) {

		d.a[0].output(d.a_ptr[shift+0], d.subset, d.n_blocks);
		d.a[1].input(d.a_ptr[shift+1], d.subset, d.n_blocks);

		d.a[0].sqrt_(d.a[1]);

		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);

	}

	void run_direct(uint shift) {
		double *a[2];
		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
		a[1] = (double *)(& (d.a_ptr[shift+1][0]));
		for(uint i=0; i<array_sub_size * simd_block_size; ++i)
				a[0][i] = sqrt(a[1][i]);
	}

};


//template <class A>
//struct Condition {
//	Expr<A> d;
//	static const uint n_flop = 3;
//
//	void run_array(uint shift) {
//
//		d.a[0].output(d.a_ptr[shift+0], d.subset, d.n_blocks);
//		d.a[1].input(d.a_ptr[shift+1], d.subset, d.n_blocks);
//		d.a[2].input(d.a_ptr[shift+2], d.subset, d.n_blocks);
//		d.a[3].input(d.a_ptr[shift+3], d.subset, d.n_blocks);
//		d.a[4].middle(d.temp[0], d.flat_subset, d.n_blocks);
//		d.a[5].middle(d.temp[1], d.flat_subset, d.n_blocks);
//		d.a[6].middle(d.temp[2], d.flat_subset, d.n_blocks);
//
//
//		d.a[0].sqrt_(d.a[1]);
//
//		d.a[0].copy_out(d.a_ptr[shift+0], d.subset, d.n_blocks);
//
//	}
//
//	void run_direct(uint shift) {
//		double *a[4];
//		a[0] = (double *)(& (d.a_ptr[shift+0][0]));
//		a[1] = (double *)(& (d.a_ptr[shift+1][0]));
//		a[2] = (double *)(& (d.a_ptr[shift+2][0]));
//		a[3] = (double *)(& (d.a_ptr[shift+3][0]));
//		for(uint i=0; i<array_sub_size * simd_block_size; ++i)
//				a[0][i] = (a[1][i] > 12) ? a[2]  - a[3] : a[2] + a[3];
//	}
//
//};


int main() {
	{
	Sum<SubsetArray> subset;
	measure(subset);
	Sum<FlatArray> flat;
	measure(flat);
	}

	{
	Add<SubsetArray> subset;
	measure(subset);
	Add<FlatArray> flat;
	measure(flat);
	}

	{
	AddScalar<SubsetArray> subset;
	measure(subset);
	AddScalar<FlatArray> flat;
	measure(flat);
	}

	{
	NormSqr<SubsetArray> subset;
	measure(subset);
	NormSqr<FlatArray> flat;
	measure(flat);
	}

	{
	Sqrt<SubsetArray> subset;
	measure(subset);
	Sqrt<FlatArray> flat;
	measure(flat);
	}

}
