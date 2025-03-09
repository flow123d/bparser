/*
 * test_simd.cc
 *
 *  Created on: Jan 3, 2020
 *      Author: jb
 */

/**
 * Test to reproduce Seg fault problems observed with SIMD.
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <malloc.h>
#include "assert.hh"
#include "config.hh"
#include "crosscompile.hh"

typedef unsigned int uint;
const uint simd_block_size=4;
const uint array_max_size = 5;
const uint array_sub_size = 3;
typedef double double4 __attribute__((__vector_size__(32)));
//typedef std::array<double, simd_block_size> double4;




int main() {
	double4 a;
	for(uint j = 0; j < simd_block_size; ++j)
		a[j] = 1;

	double4 b[2];
	for(uint i=0; i<2; ++i)
		for(uint j = 0; j < simd_block_size; ++j)
			b[i][j] = 1;

	double4 * base;
	base = b;
	for(uint i=0; i<2; ++i)
		for(uint j = 0; j < simd_block_size; ++j)
			base[i][j] = 3;
	
	base = (double4 *)align_alloc(sizeof(double) * simd_block_size * 4, sizeof(double) * simd_block_size); //no free? - LV
	base[0][0] = 1.2;

//	for(uint i=0; i<4; ++i)
//		for(uint j = 0; j < simd_block_size; ++j)
//			base[i][j] = 5;

	for(uint i=0; i<4; ++i)
		for(uint j = 0; j < simd_block_size; ++j)
			std::cout << b[i][j] << " ";


	//std::cout << "Expr constr" << "\n";
//	uint base_size = 4 * array_max_size;
//	double4 * base;
//	base = (double4 *)memalign(simd_block_size, sizeof(double) * simd_block_size * base_size);
//	ASSERT(base !=nullptr);
//	std::cout << "base: " << base << std::endl;
//	// Initialization
//	for(uint i=0; i<base_size; ++i) {
//		std::cout << "i: " << i << std::endl;
//		for(uint j = 0; j < simd_block_size; ++j) {
//			base[i][j] = i * simd_block_size + j;
//
//		}
//	}

//	double4 sum_ = {0,0,0,0};
//	double4 * a_ptr[2];
//	uint subset[array_sub_size] = {1,2,4};
//	uint flat_subset[array_sub_size] = {0,1,2};

	//std::cout << "here" << std::endl;

//	for(uint i = 0; i < 2; ++i) {
//		a_ptr[i] = base + i * array_max_size;
//	}
//
//	for(uint j=0;j<n_temp;++j) {
//		double cc = 3.14 * j;
//		for(uint k=0; k<simd_block_size; ++k)
//			constants[j][k] = cc;
//	}
//
//	for(uint j=0;j<4;++j)
//		temp[j] = base + (2 * n_arrays + j) * array_max_size;
//
//	uint seed = 1234;
//
//	make_subset(seed, subset);
//	for(uint i=0; i<n_blocks; ++i)
//		flat_subset[i] = i;
//	std::cout << "here" << std::endl;

}
