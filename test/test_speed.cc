/*
 * test_parser_speed.cc
 *
 *  Created on: Feb 2, 2020
 *      Author: jb
 */




/*
 * test_parser.cc
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */

#include <string>
#include <chrono>
#include <cmath>
#include "assert.hh"
#include "parser.hh"
#include "test_tools.hh"

#include "arena_alloc.hh"

// Optimized structure, holds data in common arena
struct ExprData {
	ExprData(uint vec_size, uint simd_size)
	: vec_size(vec_size), simd_size(simd_size)
	{
		uint simd_bytes = sizeof(double) * simd_size;

		arena = std::make_shared<bparser::ArenaAlloc>(simd_bytes, 512 * 1012);
		v1 = arena->create_array<double>(vec_size * 3);
		fill_seq(v1, 100, 100 + 3 * vec_size);
		v2 = arena->create_array<double>(vec_size * 3);
		fill_seq(v2, 200, 200 + 3 * vec_size);
		v3 = arena->create_array<double>(vec_size * 3);
		fill_seq(v3, 300, 300 + 3 * vec_size);
		v4 = arena->create_array<double>(vec_size * 3);
		fill_seq(v4, 400, 400 + 3 * vec_size);
		vres = arena->create_array<double>(vec_size * 3);
		fill_const(vres, 3 * vec_size, -100);
		subset = arena->create_array<uint>(vec_size);
		for(uint i=0; i<vec_size/simd_size; i++) subset[i] = i;
		cs1 = 4;
		for (uint i = 0; i < 3; i++)
		{
			cv1[i] = (i+1)*3;
		}
	}	

	~ExprData()
	{}

	std::shared_ptr<bparser::ArenaAlloc> arena;
	uint vec_size;
	uint simd_size;
	double *v1, *v2, *v3, *v4, *vres;
	double cs1;
	double cv1[3];
	uint *subset;

};



// Unoptimized structure, holds data in separated arrays, copies data to arenas
struct ExprData2 {
	ExprData2(uint vec_size, uint simd_size)
	: arena_1(sizeof(double)*simd_size, 512 * 1012), arena_2(sizeof(double)*simd_size, 512 * 1012), arena_3(sizeof(double)*simd_size, 512 * 1012), arena_4(sizeof(double)*simd_size, 512 * 1012),
	  arena_res(sizeof(double)*simd_size, 512 * 1012), arena_subs(sizeof(double)*simd_size, 512 * 1012), vec_size(vec_size), simd_size(simd_size)
	{
		d1 = new double[3 * vec_size];
		fill_seq(d1, 100, 100 + 3 * vec_size);
		v1 = arena_1.create_array<double>(vec_size * 3);
		d2 = new double[3 * vec_size];
		fill_seq(d2, 200, 200 + 3 * vec_size);
		v2 = arena_2.create_array<double>(vec_size * 3);
		d3 = new double[3 * vec_size];
		fill_seq(d3, 300, 300 + 3 * vec_size);
		v3 = arena_3.create_array<double>(vec_size * 3);
		d4 = new double[3 * vec_size];
		fill_seq(d4, 400, 400 + 3 * vec_size);
		v4 = arena_4.create_array<double>(vec_size * 3);
		dres = new double[3 * vec_size];
		fill_const(dres, 3 * vec_size, -100);
		vres = arena_res.create_array<double>(vec_size * 3);
		subset = arena_subs.create_array<uint>(vec_size);
		for(uint i=0; i<vec_size/simd_size; i++) subset[i] = i;
		cs1 = 4;
	}

	~ExprData2()
	{
		delete[] d1;
		delete[] d2;
		delete[] d3;
		delete[] d4;
		delete[] dres;
	}

	void copy_to_arena() {
		for (uint i=0; i<3*vec_size; ++i) {
			v1[i] = d1[i];
			v2[i] = d2[i];
			v3[i] = d3[i];
			v4[i] = d4[i];
		}
	}

	void copy_from_arena() {
		for (uint i=0; i<3*vec_size; ++i) {
			dres[i] = vres[i];
		}
	}

	bparser::ArenaAlloc arena_1;       ///< Arena of data vector 1
	bparser::ArenaAlloc arena_2;       ///< Arena of data vector 2
	bparser::ArenaAlloc arena_3;       ///< Arena of data vector 3
	bparser::ArenaAlloc arena_4;       ///< Arena of data vector 4
	bparser::ArenaAlloc arena_res;     ///< Arena of result vector
	bparser::ArenaAlloc arena_subs;    ///< Arena of subset vector
	uint vec_size;
	uint simd_size;
	double *d1, *d2, *d3, *d4, *dres;  ///< Data initialized out of arenas
	double *v1, *v2, *v3, *v4, *vres;  ///< Pointers to arenas
	double cs1;
	double cv1[3];
	uint *subset;

};




// C++ evaluation of expression "v1 + v2 + v3 + v4"
void expr1(ExprData &data) {
	for(uint i_comp=0; i_comp < 3*data.vec_size; i_comp += data.vec_size) {
		for(uint i=0; i<data.vec_size/data.simd_size; ++i) {
			uint j = i_comp + data.simd_size*data.subset[i];
			for(uint k = 0; k<data.simd_size; k++) {
				double v1 = data.v1[j+k];
				double v2 = data.v2[j+k];
				double v3 = data.v3[j+k];
				double v4 = data.v4[j+k];
				data.vres[j+k] = v1 + v2 + v3 + v4 ;
			}
		}
	}
}


// C++ evaluation of expression "3 * v1 + cs1 * v2 + v3 + 2.5 * v4"
void expr2(ExprData &data) {
	for(uint i_comp=0; i_comp < 3*data.vec_size; i_comp += data.vec_size) {
		for(uint i=0; i<data.vec_size/data.simd_size; ++i) {
			uint j = i_comp + data.simd_size*data.subset[i];
			for(uint k = 0; k<data.simd_size; k++) {
				double v1 = data.v1[j+k];
				double v2 = data.v2[j+k];
				double v3 = data.v3[j+k];
				double v4 = data.v4[j+k];
				data.vres[j+k] = 3 * v1  + data.cs1 * v2 + v3 + 2.5 * v4 ;
			}
		}
	}
}


// C++ evaluation of expression "sin(v1)"
void expr3(ExprData &data) {
	for(uint i_comp=0; i_comp < 3*data.vec_size; i_comp += data.vec_size) {
		for(uint i=0; i<data.vec_size/data.simd_size; ++i) {
			uint j = i_comp + data.simd_size*data.subset[i];
			for(uint k = 0; k<data.simd_size; k++) {
				double v1 = data.v1[j+k];
				//double v2 = data.v2[j+k];
				//double v3 = data.v3[j+k];
				//double v4 = data.v4[j+k];
				data.vres[j+k] = sin(v1);
			}
		}
	}
}


// C++ evaluation of expression "[v2, v2, v2] @ v1 + v3"
void expr4(ExprData &data) {
	for(uint i_comp=0; i_comp < 3*data.vec_size; i_comp += data.vec_size) {
		for(uint i=0; i<data.vec_size/data.simd_size; ++i) {
			uint j = i_comp + data.simd_size*data.subset[i];
			for(uint k = 0; k<data.simd_size; k++) {
				double v1 = data.v1[j+k];
				double v2 = data.v2[j+k];
				double v3 = data.v3[j+k];
				double v4 = data.v4[j+k];
				data.vres[j+k] = std::min(v1, v2) + std::max(v3, v4);
				//expression will be changed
			}
		}
	}
}


/**
 * @param expr       Parser expression
 * @param block_size Number of floats
 * @param i_expr     Specifies C++ expression function
 */
void test_expr(std::string expr, uint block_size, void (* func)(ExprData&)) {
	using namespace bparser;
	uint vec_size = 1*block_size;

	uint simd_size = get_simd_size();

	// TODO: allow changing variable pointers, between evaluations
	// e.g. p.set_variable could return pointer to that pointer
	// not so easy for vector and tensor variables, there are many pointers to set
	// Rather modify the test to fill the
	uint n_repeats = (1024 / block_size) * 100000;

	ExprData  data1(vec_size, simd_size);
	ExprData2 data2(vec_size, simd_size);
	ExprData  data3(vec_size, simd_size);
	ExprData  data4(2*vec_size, simd_size);

	double parser_time_optim, parser_time_shared_arena, parser_time_copy, parser_time_noopt, cpp_time;

	{ // one allocation in common arena
		Parser p(block_size);
		p.parse(expr);
		p.set_constant("cs1", {}, 	{data1.cs1});
		p.set_constant("cv1", {3}, 	std::vector<double>(data1.cv1, data1.cv1+3));
		p.set_variable("v1", {3}, data1.v1);
		p.set_variable("v2", {3}, data1.v2);
		p.set_variable("v3", {3}, data1.v3);
		p.set_variable("v4", {3}, data1.v4);
		p.set_variable("_result_", {3}, data1.vres);
		//std::cout << "vres: " << vres << ", " << vres + block_size << ", " << vres + 2*vec_size << "\n";
		//std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
		//std::cout.flush();
		p.compile();

		std::vector<uint> ss = std::vector<uint>(data1.subset, data1.subset+vec_size/simd_size);
		p.set_subset(ss);
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
			p.run();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		parser_time_optim = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
	}

	{ // one allocation in common arena, set this arena to processor
		Parser p(block_size);
		p.parse(expr);
		p.set_constant("cs1", {}, 	{data1.cs1});
		p.set_constant("cv1", {3}, 	std::vector<double>(data1.cv1, data1.cv1+3));
		p.set_variable("v1", {3}, data1.v1);
		p.set_variable("v2", {3}, data1.v2);
		p.set_variable("v3", {3}, data1.v3);
		p.set_variable("v4", {3}, data1.v4);
		p.set_variable("_result_", {3}, data1.vres);
		//std::cout << "vres: " << vres << ", " << vres + block_size << ", " << vres + 2*vec_size << "\n";
		//std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
		//std::cout.flush();
		p.compile(data1.arena);

		std::vector<uint> ss = std::vector<uint>(data1.subset, data1.subset+vec_size/simd_size);
		p.set_subset(ss);
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
			p.run();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		parser_time_shared_arena = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
	}

	// { // one allocation in common arena, use set_var_copy
	// 	Parser p(block_size);
	// 	p.parse(expr);
	// 	p.set_constant("cs1", {}, 	{data4.cs1});
	// 	p.set_constant("cv1", {3}, 	std::vector<double>(data4.cv1, data4.cv1+3));
	// 	p.set_var_copy("v1", {3}, data4.v1);
	// 	p.set_var_copy("v2", {3}, data4.v2);
	// 	p.set_var_copy("v3", {3}, data4.v3);
	// 	p.set_var_copy("v4", {3}, data4.v4);
	// 	p.set_variable("_result_", {3}, data4.vres);
	// 	//std::cout << "vres: " << vres << ", " << vres + block_size << ", " << vres + 2*vec_size << "\n";
	// 	//std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
	// 	//std::cout.flush();
	// 	p.compile();

	// 	std::vector<uint> ss = std::vector<uint>(data4.subset, data4.subset+vec_size/simd_size);
	// 	p.set_subset(ss);
	// 	auto start_time = std::chrono::high_resolution_clock::now();
	// 	for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
	// 		p.run();
	// 	}
	// 	auto end_time = std::chrono::high_resolution_clock::now();
	// 	parser_time_copy = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
	// }

	{ // unoptimized allocation in separated arenas
		Parser p(block_size);
		p.parse(expr);
		p.set_constant("cs1", {}, 	{data2.cs1});
		p.set_constant("cv1", {3}, 	std::vector<double>(data2.cv1, data2.cv1+3));
		p.set_variable("v1", {3}, data2.v1);
		p.set_variable("v2", {3}, data2.v2);
		p.set_variable("v3", {3}, data2.v3);
		p.set_variable("v4", {3}, data2.v4);
		p.set_variable("_result_", {3}, data2.vres);
		//std::cout << "vres: " << vres << ", " << vres + block_size << ", " << vres + 2*vec_size << "\n";
		//std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
		//std::cout.flush();
		p.compile();

		std::vector<uint> ss = std::vector<uint>(data2.subset, data2.subset+vec_size/simd_size);
		p.set_subset(ss);
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
			data2.copy_to_arena();
			p.run();
			data2.copy_from_arena();
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		parser_time_noopt = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
	}

	{ // C++ expression
		auto start_time = std::chrono::high_resolution_clock::now();
		for(uint i_rep=0; i_rep < n_repeats; i_rep++) func(data3);
		auto end_time = std::chrono::high_resolution_clock::now();
		cpp_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
	}

	// check
	double p_sum = 0;
	double c_sum = 0;
	double diff = 0;
	for(uint dim=0; dim < 3; dim++) {
		for(uint i=0; i<data1.vec_size; i++) {
			double v1 = data1.vres[dim*data1.vec_size + i];
			double v2 = data1.vres[dim*data2.vec_size + i];
			//std::cout << "res: " << v1 <<std::endl;
			diff += std::fabs(v1 - v2);
			p_sum += v1;
			c_sum += v2;
		}
	}

	std::cout << "=== Parsing of expression: '" << expr << "' ===\n";
	std::cout << "Block size: " << block_size << "\n";
	std::cout << "Diff: " << diff << " parser: " << p_sum << " c++: " << c_sum << "\n";
	std::cout << "parser time optim   : " << parser_time_optim << "\n";
	std::cout << "parser shared arena : " << parser_time_shared_arena << "\n";
	std::cout << "parser time copy    : " << parser_time_copy << "\n";
	std::cout << "parser time noopt   : " << parser_time_noopt << "\n";
	std::cout << "c++ time            : " << cpp_time << "\n";
	std::cout << "fraction: " << parser_time_optim/cpp_time << "\n";
	double n_flop = n_repeats * vec_size * 9;
	std::cout << "parser FLOPS: " << n_flop / parser_time_optim << "\n";
	std::cout << "c++ FLOPS   : " << n_flop / cpp_time << "\n";
	std::cout << "======================================================\n\n";

}




void test_expression() {
	std::vector<uint> block_sizes = {64, 256, 1024};
	for (uint i=0; i<block_sizes.size(); ++i) {
		test_expr("v1 + v2 + v3 + v4", block_sizes[i], &expr1);
		test_expr("3 * v1 + cs1 * v2 + v3 + 2.5 * v4", block_sizes[i], &expr2);
		test_expr("sin(v1)", block_sizes[i], &expr3);
		test_expr("[v2, v2, v2] @ v1 + v3", block_sizes[i], &expr4); // correct expression
	}
}



int main()
{
	test_expression();
}




