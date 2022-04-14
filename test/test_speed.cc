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
#include "assert.hh"
#include "parser.hh"
#include "test_tools.hh"
#include "processor.hh"
#include "vectorclass.h"

#include "arena_alloc.hh"



using namespace bparser;


uint get_simd_size()
{
	if (__builtin_cpu_supports("avx512f"))
	{
		return 8;
	}
	if (__builtin_cpu_supports("avx2"))
	{
		return 4;
	}
	if (__builtin_cpu_supports("avx"))
	{
		return 4;
	}
	if (__builtin_cpu_supports("sse"))
	{
		return 2;
	}
	else
	{
		return 1;
	}
}


struct ExprData {
	ExprData(uint vec_size, uint simd_size)
	: vec_size(vec_size)
	{
		ExprData::simd_size = simd_size;
		uint simd_bytes = sizeof(double) * simd_size;
		arena = new ArenaAlloc(simd_bytes, var_cnt * vec_size * sizeof(double));


		// std::cout << "\nIn test_speed.cc, ExprData:" << std::endl;
		// std::cout << "var_cnt: " << var_cnt << std::endl;
		// std::cout << "arena: " << var_cnt * vec_size * sizeof(double) << std::endl;
		
		
		v1 = (*arena).create_array<double>(vec_size * 3);
		fill_seq(v1, 100, 100 + 3 * vec_size);

		v2 = (*arena).create_array<double>(vec_size * 3);
		fill_seq(v2, 200, 200 + 3 * vec_size);
		
		v3 = (*arena).create_array<double>(vec_size * 3);
		fill_seq(v3, 100, 100 + 3 * vec_size * 0.5, 0.5);
		
		v4 = (*arena).create_array<double>(vec_size * 3);
		fill_seq(v4, 200, 200 + 3 * vec_size * 0.5, 0.5);
		
		vres = (*arena).create_array<double>(vec_size * 3);
		fill_const(vres, 3 * vec_size, -100);
		
		subset = (*arena).create_array<uint>(vec_size);
		for(uint i=0; i<vec_size/simd_size; i++) subset[i] = i;

		cs1 = 4;

		for (uint i = 0; i < 3; i++)
		{
			cv1[i] = (i+1)*3;
		}
	}

	~ExprData()
	{
		delete arena;
	}

	uint vec_size;
	uint simd_size;
	double *v1, *v2, *v3, *v4, *vres;
	double cs1;
	double cv1[3];
	uint *subset;
	ArenaAlloc *arena;


	uint vec_cnt = 4;
	uint cvec_cnt = 1;
	uint const_cnt = 1;

	uint var_cnt = 3 * (vec_cnt + cvec_cnt + 1) + const_cnt + 1;  	// +1 je na result (vres) a na subset
};



void expr1(ExprData &data) {
	for(uint i_comp=0; i_comp < 3*data.vec_size; i_comp += data.vec_size) {
		for(uint i=0; i<data.vec_size/data.simd_size; ++i) {
			uint j = i_comp + data.simd_size*data.subset[i];
			for(uint k = 0; k<data.simd_size; k++) {
				double v1 = data.v1[j+k];
				double v2 = data.v2[j+k];

				// data.vres[j+k] = v1 * v2;
				// data.vres[j+k] = v1 * v2 * v1;
				data.vres[j+k] = v1 * v2 * v1 * v2;

				// double v3 = data.v3[j+k];
				// data.vres[j+k] = v1 * v2 * v3;

				// double v3 = data.v3[j+k];
				// double v4 = data.v4[j+k];
				// double cv1 = data.cv1[j+k];
				// data.vres[j+k] = 3 * v1 + data.cs1 * v2 + cv1 * v3 + pow(v4,2);
			}
		}
	}
}


void test_expr(std::string expr) {
	using namespace bparser;
	uint block_size = 1024; // number of floats
	uint vec_size = 1*block_size;

	// TODO: allow changing variable pointers, between evaluations
	// e.g. p.set_variable could return pointer to that pointer
	// not so easy for vector and tensor variables, there are many pointers to set
	// Rather modify the test to fill the
	uint n_repeats = 1000;
	uint simd_size = get_simd_size();

	
	ExprData data1(vec_size, simd_size);
	ExprData data2(vec_size, simd_size);


	Parser p(block_size);
	p.parse(expr);
	p.set_constant("cs1", {}, {data1.cs1});
	p.set_constant("cv1", {3}, std::vector<double>(data1.cv1, data1.cv1+3));
	p.set_variable("v1", {3}, data1.v1);
	p.set_variable("v2", {3}, data1.v2);
	p.set_variable("v3", {3}, data1.v3);
	p.set_variable("v4", {3}, data1.v4);
	p.set_variable("_result_", {3}, data1.vres);

	//std::cout << "vres: " << vres << ", " << vres + block_size << ", " << vres + 2*vec_size << "\n";
	//std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
	//std::cout.flush();
	ExpressionDAG se = p.compile();

	ProcessorBase * processor = create_processor(se, vec_size, simd_size);
	p.set_processor(processor);

	std::vector<uint> ss = std::vector<uint>(data1.subset, data1.subset + vec_size / simd_size); //bylo lomeno 4
	p.set_subset(ss);
	auto start_time = std::chrono::high_resolution_clock::now();
	for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
		p.run();
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	double parser_time  =
			std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

	start_time = std::chrono::high_resolution_clock::now();
	for(uint i_rep=0; i_rep < n_repeats; i_rep++) {
		expr1(data2);
	}
	end_time = std::chrono::high_resolution_clock::now();
	double cpp_time  =
			std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

	// check
	double p_sum = 0;
	double c_sum = 0;
	double diff = 0;
	for(uint dim=0; dim < 3; dim++) {
		for(uint i=0; i<data1.vec_size; i++) {
			double v1 = data1.vres[dim*data1.vec_size + i];
			double v2 = data2.vres[dim*data2.vec_size + i];
			//std::cout << "res: " << v1 <<std::endl;
			diff += std::fabs(v1 - v2);
			p_sum += v1;
			c_sum += v2;
		}
	}
	arena_1.destroy();
	arena_2.destroy();
	std::cout << "In test_speed.cc" << std::endl;
	std::cout << "Diff: " << diff << " parser: " << p_sum << " c++: " << c_sum << "\n";
	std::cout << "parser time : " << parser_time << "\n";
	std::cout << "c++ time    : " << cpp_time << "\n";
	std::cout << "fraction: " << parser_time/cpp_time << "\n";
	double n_flop = n_repeats * vec_size * 9;
	std::cout << "parser FLOPS: " << n_flop / parser_time << "\n";
	std::cout << "c++ FLOPS   : " << n_flop / cpp_time << "\n";
	
	// std::cout << "velikost Vec2d        : " << sizeof(Vec2d) << "\n";
	// std::cout << "velikost Vec4d        : " << sizeof(Vec4d) << "\n";
	// std::cout << "velikost Vec8d        : " << sizeof(Vec8d) << "\n";
	// std::cout << "velikost double       : " << sizeof(double) << "\n";
	// std::cout << "velikost Proc<Vec2d>  : " << sizeof(Processor<Vec<Vec2d>>) << "\n";
	// std::cout << "velikost Proc<Vec4d>  : " << sizeof(Processor<Vec<Vec4d>>) << "\n";
	// std::cout << "velikost Proc<Vec8d>  : " << sizeof(Processor<Vec<Vec8d>>) << "\n";
	// std::cout << "velikost Proc<double> : " << sizeof(Processor<Vec<double>>) << "\n";
}




void test_expression() {
	// test_expr("3 * v1 + cs1 * v2 + cv1 * v3 + v4**2");
	// test_expr("v1 * v2");
	// test_expr("v1 * v2 * v1");
	// test_expr("v1 * v2 * v3");
	test_expr("v1 * v2 * v1 * v2");
}



int main()
{
	test_expression();
}




