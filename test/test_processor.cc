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
	double v2_values[vec_size]; // 5 * double 4
	double v3_values[vec_size]; // 5 * double 4

	ScalarExpression se;
	ScalarNode * c1 = se.create_const(3.14);
	ScalarNode * v1 = se.create_value(v1_values);
	ScalarNode * e1 = se.create<_abs_>(v1);
	ScalarNode * e2 = se.create<_add_>(e1, c1);
	ScalarNode * r1 = se.create_result(e1, v2_values);
	ScalarNode * r2 = se.create_result(e2, v3_values);

	std::vector<ScalarNode *> results = {r1, r2};
	Processor * processor = Processor::create(se, vec_size);


	std::vector<uint> subset = {1, 3, 4};
	processor->set_subset(subset);
	processor->run();


	subset = {0, 2, 4};
	processor->set_subset(subset);
	processor->run();

}

