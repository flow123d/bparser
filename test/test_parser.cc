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
	auto vars = p.free_symbols();

	std::cout << "free vars test: " << expr << "\n";
	bool success = (vars == ref_vars);
	if (!success) {
		std::cout <<  "vars:  ";
		for(std::string v : vars) std::cout << v << ", ";
		std::cout <<  "ref_vars:  ";
		for(std::string v : ref_vars) std::cout << v << ", ";
		std::cout << "\n";
		std::cout.flush();
	}
	return success;
}


void test_free_variables() {
	std::cout << "\n" << "** test free variables" << "\n";
	EXPECT(test_fv("1+2", {}));
	EXPECT(test_fv("pi+e+epsilon+dd", {"dd"}));
	EXPECT(test_fv("a+b", {"a", "b"}));
	EXPECT(test_fv("a=1;a+b", {"b"}));
}

constexpr uint vec_size = 8;

std::vector<double> eval_expr_(std::string expr, bparser::Shape ref_shape = {}) {
	std::cout << "parser test : " << expr << "\n";
	using namespace bparser;

	//auto m1 = new double[vec_size * 6];
	std::vector<double> as1(vec_size, 1);
	std::vector<double> av2(3*vec_size, 2);

	Parser p(vec_size);
	p.parse(expr);
	std::cout << "  AST: " << p.print_ast() << "\n";
	p.set_variable("as1", {}, &(as1[0]));
	p.set_variable("av2", {3}, &(av2[0]));
        
//        auto mat = new double[vec_size * 3 * 3]; // (0,0,0),...,(0,0,999), (0,1,0),...,(0,1,999), ...
//        p.set_variable("mat", {3, 3}, mat);
        
	p.set_constant("cs3", {}, 	{3});
	p.set_constant("cv4", {3}, 	{4, 5, 6});
	//p.set_variable("_result_", {result_size}, vres);
//	std::cout << "vres: " << vres << ", " << vres + vec_size << ", " << vres + 2*vec_size << "\n";
//	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
//	std::cout.flush();
	p.compile();
	double * vres = p.tmp_result_ptr();
	Shape res_shape = p.result_array().shape();
	if (ref_shape.size() > 0 and ! bparser::same_shape(ref_shape, res_shape)) {
		std::cout << "Wrong result shape.\n" <<
				"ref: " << print_vector(ref_shape) << "\n"
				"res: " << print_vector(res_shape) << "\n";
	}

	uint result_size = shape_size(p.result_array().shape());
	fill_const(vres, vec_size * result_size, -1e100); // undefined value
	p.set_subset({0, 1});
	p.run();

	std::vector<double> res(result_size * vec_size);
	for(uint i=0; i < res.size(); i++) res[i] = vres[i];
	p.destroy_processor();
	return res;
}

bool test_expr(std::string expr, std::vector<double> ref_result, bparser::Shape ref_shape = {}) {
	auto res = eval_expr_(expr, ref_shape);

	// check
	bool success = true;
	for(uint i=0; i < vec_size; i++) {
		for(uint j=0; j < ref_result.size(); j++) {
			if (fabs(res[j * vec_size + i] - ref_result[j]) > 1e-3 * fabs(ref_result[j]) + 1e-15 ) {
				success = false;
				std::cout << "  " << i << "," << j <<
				" ref: " << ref_result[j] <<
				" res: " << res[j * vec_size + i] << "\n";
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

std::vector<double> eval_bool_expr_(std::string expr) {
	std::cout << "parser test : " << expr << "\n";
	using namespace bparser;

	uint simd_size = 1;
	uint simd_bytes = sizeof(double) * simd_size;
	ArenaAlloc *arena = new ArenaAlloc(simd_bytes, 6 * vec_size * sizeof(double));

	auto v1 = (*arena).create_array<double>(vec_size * 3);
	fill_seq(v1, 88, 100 + 3 * vec_size * 2, 2.0);

	auto v2 = (*arena).create_array<double>(vec_size * 3);
	fill_seq(v2, 100, 100 + 3 * vec_size);

	Parser p(vec_size);
	p.parse(expr);
	std::cout << "  AST: " << p.print_ast() << "\n";
	p.set_variable("v1", {3}, v1);
	p.set_variable("v2", {3}, v2);
	// p.set_constant("cs3", {}, 	{3});
	// p.set_constant("cv4", {3}, 	{4, 5, 6});
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

bool test_bool_expr(std::string expr, std::vector<double> ref_result) {
	auto res = eval_bool_expr_(expr);

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

void test_bool_expression()
{
	std::cout << std::endl << "** test bool expression" << std::endl;
	BP_ASSERT(test_expr("v1 < v2", {1,1,0}));
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
	BP_ASSERT(test_expr("cv4[[0]] ** 2", {16}));
	//BP_ASSERT(test_expr("m=[cv4, av2]; m[[0,0,1], [0, 2, 0]]", {4, 6, 2}, {3}));

	BP_ASSERT(fail_expr("cs3[0]", "Too many indices")); // ?fail
	BP_ASSERT(test_expr("cv4[0, None] * cv4[None, 1]", {})); // ?? matrix
	//BP_ASSERT(test_expr("[]", {0}));
	BP_ASSERT(fail_expr("[]", "Empty Array"));
	BP_ASSERT(fail_expr("[1,cv4,av2]", "stack: all input arrays must have the same shape"));
	BP_ASSERT(test_expr("[1,1]", {1, 1}));
	BP_ASSERT(test_expr("a=[[1,1,1], cv4, av2]; a[:, 0]", {1, 4, 2}));
	BP_ASSERT(test_expr("[[1,2], [3,4]]", {1, 2, 3, 4}, {2,2}));

	//BP_ASSERT(test_expr("cv4[2] ** av2", {25, 25, 25}));
	//BP_ASSERT(test_expr("cv4[2] ** av2", {25, 25, 25}));

	auto vec_false = std::vector<double>(1, bparser::details::double_false());
	auto vec_true = std::vector<double>(1, bparser::details::double_true());
	BP_ASSERT(test_expr("cs3 > 4.5", vec_false));
	BP_ASSERT(test_expr("cs3 < 4.5", vec_true));

	BP_ASSERT(test_expr("2 < cs3 < 4.5", vec_true));
	BP_ASSERT(test_expr("(2 < cs3) < 2", vec_true));

	BP_ASSERT(test_expr("3 if cs3 < 4.5 else 4", {3}));
	BP_ASSERT(test_expr("3 if cs3 > 4.5 else 4", {4}));


	BP_ASSERT(test_expr("[3, 4] @ [[1], [2]]", {11}));
	BP_ASSERT(test_expr("[3, 4, 1] @ [[1], [2], [3]]", {14}));
	BP_ASSERT(test_expr("[[1, 2], [2, 3], [3, 4]] @ [[1], [2]]", {5, 8, 11}));
	//BP_ASSERT(test_expr("0 if cv4 > 4.5 else 1", {1, 0, 0}));

	BP_ASSERT(test_expr("minimum([1,2,3], [0,4,3])", {0,2,3}));
	BP_ASSERT(test_expr("maximum([1,2,3], [0,4,3])", {1,4,3}));

	BP_ASSERT(test_expr("sin(pi*as1)", {0})); // sin(pi*1)

	BP_ASSERT(test_expr("flatten([[1,2],[3,4]])", {1, 2, 3, 4}, {4}));
	BP_ASSERT(test_expr("eye(2)[0,0]", {1}));
	BP_ASSERT(test_expr("eye(2)[1,0]", {0}));
	BP_ASSERT(test_expr("eye(2)", {1, 0, 0, 1}, {2, 2}));
	BP_ASSERT(test_expr("flatten(eye(3))", {1, 0, 0, 0, 1, 0, 0, 0, 1}));
	BP_ASSERT(test_expr("zeros([2,3])", {0, 0, 0, 0, 0, 0}, {2,3}));
	BP_ASSERT(test_expr("ones([2,3])", {1, 1, 1, 1, 1, 1}, {2,3}));
	//BP_ASSERT(test_expr("full([2,3], 5)", {5, 5, 5, 5, 5, 5}, {2,3}));
	//BP_ASSERT(test_expr("norm([2, 3])", {5}));
}


void test_speed_cases() {

}

int main()
{
	// test_free_variables();
	test_expression();
	// test_bool_expression();
#ifdef NDEBUG
	test_speed_cases();
#endif
}




