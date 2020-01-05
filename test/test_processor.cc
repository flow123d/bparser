/*
 * test_processor.cc
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */



#include <string>
#include "assert.hh"
#include "scalar_expr.hh"
#include "processor.hh"









int main()
{
	using namespace bparser;
	using namespace bparser::details;

	const uint vec_size = 20;
	double v1_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v1_values[i] = (double)i - 10;
	//std::cout << "v1_ptr: " << v1_values << "\n";

	double v2_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v2_values[i] = -100;
	//std::cout << "v2_ptr: " << v2_values << "\n";

	double v3_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v3_values[i] = -100;
	//std::cout << "v3_ptr: " << v3_values << "\n";


	ScalarExpression se;
	ScalarNode * c1 = se.create_const(100);
	ScalarNode * v1 = se.create_value(v1_values);
	ScalarNode * e1 = se.create<_abs_>(v1);
	ScalarNode * e2 = se.create<_add_>(e1, c1);
	ScalarNode * r1 = se.create_result(e1, v2_values);
	ScalarNode * r2 = se.create_result(e2, v3_values);

	auto sorted = se.sort_nodes();
//	for(auto node : sorted) {
//		std::cout << "node: " << node << " res: " << node->result_storage << " sidx: " << node->result_idx_ << "\n";
//	}

	ASSERT(sorted[5] == v1);
	ASSERT(v1->result_idx_ == 0);
	ASSERT(sorted[4] == e1);
	ASSERT(e1->result_idx_ == 1);
	ASSERT(sorted[3] == r1);
	ASSERT(r1->result_idx_ == 1);
	ASSERT(sorted[2] == c1);
	ASSERT(c1->result_idx_ == 2);
	ASSERT(sorted[1] == e2);
	ASSERT(e2->result_idx_ == 3);
	ASSERT(sorted[0] == r2);
	ASSERT(r2->result_idx_ == 3);
	Processor * processor = Processor::create(se, vec_size);


	std::vector<uint> subset = {1, 3, 4};
	processor->set_subset(subset);
	processor->run();
	// v1 active: [ -6, -5, -4, -3; 2, 3, 4, 5; 6, 7, 8, 9 ]
	// r1, e1     [  6   5   4   3  2  3  4  5  6  7  8  9 ]
	// r2, e2     [  6   5   4   3  2  3  4  5  6  7  8  9 ]  + 100
	for(uint i=0;i < vec_size; ++i)
		std::cout << "i: " << i << " v2: " << v2_values[i]  << " v3: " <<v3_values[i] << "\n";

	ASSERT(v2_values[0] == -100);
	ASSERT(v2_values[3] == -100);
	ASSERT(v2_values[4] ==  6);
	ASSERT(v2_values[7] ==  3);
	ASSERT(v2_values[8] == -100);
	ASSERT(v2_values[11] == -100);
	ASSERT(v2_values[12] == 2);
	ASSERT(v2_values[15] == 5);
	ASSERT(v2_values[16] == 6);
	ASSERT(v2_values[19] == 9);

	ASSERT(v3_values[0] == -100);
	ASSERT(v3_values[3] == -100);
	ASSERT(v3_values[4] ==  106);
	ASSERT(v3_values[7] ==  103);
	ASSERT(v3_values[8] == -100);
	ASSERT(v3_values[11] == -100);
	ASSERT(v3_values[12] == 102);
	ASSERT(v3_values[15] == 105);
	ASSERT(v3_values[16] == 106);
	ASSERT(v3_values[19] == 109);

//	subset = {0, 2, 4};
//	processor->set_subset(subset);
//	processor->run();

}

