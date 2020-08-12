/*
 * test_parser.cc
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */

#include <string>

#include "test_tools.hh"
#include "assert.hh"
#include "parser.hh"



bool test_fv(std::string expr, std::vector<std::string> ref_vars) {
	using namespace bparser;
	Parser p(4);
	p.parse(expr);
	auto vars = p.variables();

	std::cout << "free vars test: " << expr << "\n";
	bool success = (vars == ref_vars);
	if (!success) {
		std::cout <<  "  ";
		for(std::string v : vars)
			std::cout << v << ", ";
		std::cout << "\n";
		std::cout.flush();
	}
	return success;
}


void test_free_variables() {
	std::cout << "\n" << "** test free variables" << "\n";
	EXPECT(test_fv("1+2", {}));
	EXPECT(test_fv("a+b", {"a", "b"}));
	EXPECT(test_fv("a=1;a+b", {"b"}));
}

constexpr uint vec_size = 8;

std::vector<double> eval_expr_(std::string expr) {
	std::cout << "parser test : " << expr << "\n";
	using namespace bparser;

	//auto m1 = new double[vec_size * 6];
	auto as1 = new double[vec_size];
	fill_const(as1, vec_size, 1);
	auto av2 = new double[vec_size * 3];
	fill_const(av2, 3 * vec_size, 2);





	Parser p(vec_size);
	p.parse(expr);
	std::cout << "  AST: " << p.print_ast() << "\n";
	p.set_variable("as1", {}, as1);
	p.set_variable("av2", {3}, av2);
	p.set_constant("cs3", {}, 	{3});
	p.set_constant("cv4", {3}, 	{4, 5, 6});
	//p.set_variable("_result_", {result_size}, vres);
//	std::cout << "vres: " << vres << ", " << vres + vec_size << ", " << vres + 2*vec_size << "\n";
//	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
//	std::cout.flush();
	p.compile();
	double * vres = p.tmp_result_ptr();
	uint result_size = shape_size(p.result_array().shape());
	fill_const(vres, vec_size * result_size, -1e100); // undefined value
	p.set_subset({0, 1});
	p.run();

	std::vector<double> res(result_size * vec_size);
	for(uint i=0; i < res.size(); i++) res[i] = vres[i];
	return res;
}

bool test_expr(std::string expr, std::vector<double> ref_result) {
	auto res = eval_expr_(expr);

	// check
	bool success = true;
	for(uint i=0; i < vec_size; i++) {
		for(uint j=0; j < ref_result.size(); j++) {
			if (fabs(res[j * vec_size + i] - ref_result[j]) > 1e-3 * fabs(ref_result[j]) ) {
				success = false;
				std::cout << "  " << i << "," << j <<
				" ref: " << ref_result[j] <<
				" res: " << res[i * ref_result.size() + j] << "\n";
			}
		}
	}
	std::cout.flush();
	return success;
}

bool fail_expr(std::string expr, std::string ref_msg) {
	try {
		eval_expr_(expr);
	} catch (bparser::Exception &e) {
		size_t pos = std::string(e.what()).find(ref_msg);
		if (pos == std::string::npos) {
			std::cout << "\nBP exception msg: " << e.what() << "\n";
			std::cout << "Expected msg: " << ref_msg << "\n";
			return false;
		}
		return true;
	}
	return false;
}


void test_expression() {
	/**
	 * All tests have defined:
	 * as1 - scalar array == 1
	 * av2 - vector array == [2, 2, 2]
	 * cs3 - scalar constant == 3
	 * cv4 - vector constant = [4,5,6]
	 */
	std::cout << "\n" << "** test expression" << "\n";
	BP_ASSERT(test_expr("cs3 * av2", {6,6,6}));
	BP_ASSERT(test_expr("cv4 * av2", {8,10,12}));
	BP_ASSERT(test_expr("as1 - cv4", {-3,-4,-5}));
	BP_ASSERT(test_expr("[2,3,4] / av2", {1,1.5,2}));


	BP_ASSERT(test_expr("cv4[1] ** av2", {25, 25, 25}));
	BP_ASSERT(test_expr("cv4[:2] ** 2", {16, 25}));
	BP_ASSERT(test_expr("cv4[[0,1]] ** 2", {16, 25}));
	BP_ASSERT(fail_expr("cs3[0]", "Too many indices")); // ?fail
	BP_ASSERT(test_expr("cv4[0, None] * cv4[None, 1]", {})); // ?? matrix
	BP_ASSERT(fail_expr("[]", "stack: need at least one array"));
	BP_ASSERT(fail_expr("[1,cv4,av2]", "stack: all input arrays must have the same shape"));
	BP_ASSERT(test_expr("[1,1]", {1, 1}));
	BP_ASSERT(test_expr("a=[[1,1,1], cv4, av2]; a[:, 0]", {1, 4, 2}));
	//BP_ASSERT(test_expr("cv4[2] ** av2", {25, 25, 25}));
	//BP_ASSERT(test_expr("cv4[2] ** av2", {25, 25, 25}));

	auto vec_false = std::vector<double>(1, bparser::details::double_false);
	auto vec_true = std::vector<double>(1, bparser::details::double_true);
	BP_ASSERT(test_expr("cs3 > 4.5", vec_false));
	BP_ASSERT(test_expr("cs3 < 4.5", vec_true));

	BP_ASSERT(test_expr("3 if cs3 < 4.5 else 4", {3}));
	BP_ASSERT(test_expr("3 if cs3 > 4.5 else 4", {4}));

	//BP_ASSERT(test_expr("cv4 > 4.5", {0, 1, 1}));
	//BP_ASSERT(test_expr("0 if cv4 > 4.5 else 1", {1, 0, 0}));
}


void test_speed_cases() {

}

int main()
{
	//test_free_variables();
	test_expression();
#ifdef NDEBUG
	test_speed_cases();
#endif
}




