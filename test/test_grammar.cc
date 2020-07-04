/*
 * test_parser.hh
 *
 *  Created on: Dec 26, 2019
 *      Author: jb
 */

#ifndef TEST_TEST_PARSER_CC_
#define TEST_TEST_PARSER_CC_

#include <string>
#include "grammar.hh"
#include "ast.hh"

#include "assert.hh"
#include "test_tools.hh"





/**
 * Return true if the string match the grammar.
 */
bool match(std::string s, std::string ref_ast) {
	std::cout << "*";
	try {
		bparser::ast::operand ast;
		bparser::parse_expr(s, ast);

		// print AST
		std::string s =  bparser::ast::print(ast);
		if (s != ref_ast) {
			std::cout << "\nParsed   AST: " << s << "\n";
			std::cout << "Expected AST: " << ref_ast << "\n";
			return false;
		}
	} catch (bparser::Exception &e) {
		std::cout << e.what() << std::endl;
		return false;
	}
	return true;
}

bool fail(std::string s, std::string ref_msg) {
	std::cout << "*";
	try {
		bparser::ast::operand ast;
		bparser::parse_expr(s, ast);

		// print AST
		std::string s =  bparser::ast::print(ast);
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

void test_primary() {
	std::cout << "\ntest_primary" << "\n";
	EXPECT(match("123", "123"));
	EXPECT(match("123.0", "123"));
	EXPECT(match("1.23e2", "123"));
	EXPECT(match("1.23e-2", "0.0123" ));
	EXPECT(match("-1.23e-2","-(0.0123)" ));
	EXPECT(match("(-1.23e-2)","-(0.0123)" ));

	EXPECT(match("sin(1.0)", "sin(1)"));
	EXPECT(fail("sin (2)", "Expected \"(\" at \" (2)\""));
	EXPECT(match("atan2(1, 2)", "atan2(1_2)"));
	EXPECT(match("e", "`e`"));

	EXPECT(match("__anything", "`__anything`"));
	EXPECT(fail("1anything", "Parsing failed at: anything"));
	EXPECT(fail("1 1", "Parsing failed at: 1"));
	EXPECT(fail("a n y", "Parsing failed at: n y"));
}

void test_operators() {
	std::cout << "\ntest_operators" << "\n";
	EXPECT(match("2 * 5", "*(2_5)"));
	EXPECT(match("-2 * 3 / 4", "/(*(-(2)_3)_4)"));
	EXPECT(match("1 * 5 / 4 + 3 - 4", "-(+(/(*(1_5)_4)_3)_4)"));
	EXPECT(match("4 + 2 * (3 - 4)", "+(4_*(2_-(3_4)))"));
	EXPECT(match("((123))", "123"));
	EXPECT(match("((123)*1)/4", "/(*(123_1)_4)"));
	EXPECT(match("(1*2) / (3*4)", "/(*(1_2)_*(3_4))"));

	EXPECT(match("2 * -5", "*(2_-(5))")); // should possibly fail
	//ASSERT(! match("-1.23e-2 * -5.3"));
	//ASSERT(! match("+ -1.23e-2"));
	EXPECT(match("+ + -1.23e-2", "+(+(-(0.0123)))"));
}

void test_arrays() {
	std::cout << "\ntest_arrays" << "\n";
	EXPECT(match("[1, 2, 3]", ",(,(,(None(0)_1)_2)_3)"));
	EXPECT(match("[[1+1, 2*1], [3/1, 4%2]]", ",(,(None(0)_,(,(None(0)_+(1_1))_*(2_1)))_,(,(None(0)_/(3_1))_%(4_2)))"));
	EXPECT(match("[[1+1, 2*1], [3/1]]", ",(,(None(0)_,(,(None(0)_+(1_1))_*(2_1)))_,(None(0)_/(3_1)))")); // should fail later in conversion of AST to Arrays
	//EXPECT(match("[1, 2, 3] + sin([xyz])", ""));
	EXPECT(fail("[(]]", ""));
	//EXPECT(fail("(1)[1]", ""));
	//EXPECT(fail("xy[1]", ""));
	//EXPECT(fail("[1,2][1]", ""));
	//EXPECT(fail("[1,x[0]]", ""));
	//EXPECT(match("[1,2][0]", ""));

}


int main()
{
	test_primary();
	test_operators();
	test_arrays();
}


#endif /* TEST_TEST_PARSER_CC_ */
