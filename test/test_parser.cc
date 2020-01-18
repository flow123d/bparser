/*
 * test_parser.cc
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */

#include <string>
#include "assert.hh"
#include "parser.hh"
#include "test_tools.hh"

void test_ast(std::string expr, std::string ref_ast) {
	using namespace bparser;
	Parser p(4);
	p.parse(expr);
	std::string s = p.print_ast();
	std::cout << "\n";
	std::cout << "Expr: " << expr << "\n";
	std::cout << "AST: " << s << "\n";
	ASSERT(s == ref_ast);
	std::cout.flush();
}

void test_ast_cases() {

	// primary rule
	// -----------

	// double
	test_ast("123.0", "123");
	// variable
	test_ast("xyz", "<xyz>");
	// constant
	test_ast("pi", "3.14159");
	// unary op
	test_ast("+sin(1)", "+(sin(1))");
	test_ast("-sin(1)", "-(sin(1))");
	// unary_fn
	test_ast("sin(123.0)", "sin(123)");
	// binary_fn
	test_ast("pow(1.2, 3.4)", "pow(1.2, 3.4)");
	// parenthesis
	test_ast("1.2 ** (3.4 ** 5.6)", "**(1.2, **(3.4, 5.6))"); // explicit
	test_ast("(1.2 ** 3.4) ** 5.6", "**(**(1.2, 3.4), 5.6)"); // explicit

	// factor rule - right associativity of power
	// -----------

	test_ast("1.2 ** 3.4", "**(1.2, 3.4)");
	test_ast("1.2 ** 3.4 ** 5.6", "**(1.2, **(3.4, 5.6))"); // yet wrong power associativity

	// multiplicative, additive, relational, equality, logical
	// -----------
	// multiplicative
	test_ast("1 * 2**1 / 3", "/(*(1, **(2, 1)), 3)");
	// additive
	test_ast("1 + 2*1 - 3", "-(+(1, *(2, 1)), 3)");
	// relational
	test_ast("1 + 2 < 2", "<(+(1, 2), 2)");
	// equality
	test_ast("1 + 2 == 2", "==(+(1, 2), 2)");
	//test_ast("1 == 2 != 3 < 4", "");
	//ASSERT_THROW(test_ast("1 == 2 != 3", ""), "Parsing failed at != 3");
	// TODO: capture all parser errors and throuw bparser exception
	// current error handler doesn't capture:
	// terminate called after throwing an instance of 'std::runtime_error'
	// what():  Parsing failed at != 3

//	// logical
//	test_ast("1 * 2**1 * 3", "bin(bin(1, bin(2, 1)), 3)");

	// program and assignment
	test_ast("a=1;a+4", ";(a = 1, +(<a>, 4))");

	// TODO: test error detection and reporting
}

void test_fv(std::string expr, std::vector<std::string> ref_vars) {
	using namespace bparser;
	Parser p(4);
	p.parse(expr);
	auto vars = p.variables();
	std::cout << "\n";
	std::cout << "Expr: " << expr << "\n";
	for(std::string v : vars)
		std::cout << v << ", ";
	std::cout << "\n";
	ASSERT(vars == ref_vars);
	std::cout.flush();

}

void test_free_variables() {
	test_fv("1+2", {});
	test_fv("a+b", {"a", "b"});
	test_fv("a=1;a+b", {"b"});
}


void test_expr(std::string expr) {
	using namespace bparser;
	uint vec_size = 8;
	double m1[vec_size * 6];
	double v1[vec_size * 3];
	fill_const(v1, 3 * vec_size, 100);
	double v2[vec_size * 3];
	fill_const(v2, 3 * vec_size, 200);
	double vres[vec_size * 3];
	fill_const(vres, 3 * vec_size, -100);

	Parser p(vec_size);
	p.parse(expr);
	p.set_constant("cs1", {}, 	{2});
	p.set_constant("cv1", {3}, 	{1, 2, 3});
	p.set_variable("v1", {3}, v1);
	p.set_variable("v2", {3}, v2);
	p.set_variable("_result_", {3}, vres);
	std::cout << "vres: " << vres << ", " << vres + vec_size << ", " << vres + 2*vec_size << "\n";
	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
	std::cout.flush();
	p.compile();
	p.set_subset({0, 1});
	p.run();
	std::cout << print_vec(vres, 3*vec_size);
}

void test_expression() {
	test_expr("1 * v1 + cs1 * v2");
}


void test_speed_cases() {

}

int main()
{
	//test_ast_cases();
	//test_free_variables();
	test_expression();
#ifdef NDEBUG
	test_speed_cases();
#endif
}




