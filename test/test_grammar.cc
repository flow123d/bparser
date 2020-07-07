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
	// literal_double
	EXPECT(match("123", "123"));
	EXPECT(match("123.0", "123"));
	EXPECT(match("1.23e2", "123"));
	EXPECT(match("1.23e-2", "0.0123" ));

	// array constructor
	// in separate test

	// enclosure
	EXPECT(match("(-1.23e-2)", "-(0.0123)" ));
	EXPECT(match("(((-1.23e-2)))", "-(0.0123)" ));
	EXPECT(match("((1 + (2 + 3)))", "+(1,+(2,3))" ));
	EXPECT(match("((1 + (2 + 3)))", "+(1,+(2,3))" ));
	EXPECT(fail("(1 + 1))", "")); // unbalanced
	EXPECT(fail("((1 + 1)", "")); // unbalanced
	EXPECT(fail("(1)(2)", "")); // missing operator

	// call
	EXPECT(match("sin(1.0)", "sin(1)"));
	EXPECT(fail("sin (2)", "Expected \"(\" at \" (2)\""));
	EXPECT(match("atan2(1, 2)", "atan2(1,2)"));

	// identifier
	EXPECT(match("e", "`e`"));
	EXPECT(match("__anything1_2", "`__anything1_2`"));
	EXPECT(fail("1anything", "Parsing failed at: anything"));
	EXPECT(fail("1 1", "Parsing failed at: 1"));
	EXPECT(fail("a n y", "Parsing failed at: n y"));

	// const_lit
	EXPECT(match("True", "True(0)"));
	EXPECT(match("False", "False(0)"));
	EXPECT(match("None", "None(0)"));

	// subscription in separate test
}


void test_arrays() {
	std::cout << "\ntest_arrays" << "\n";
	EXPECT(match("[1, 2, 3]", "array(1,2,3)"));
	EXPECT(match("[[1+1, 2*1], [3/1, 4%2]]", "array(array(+(1,1),*(2,1)),array(/(3,1),%(4,2)))"));
	EXPECT(match("[[1+1, 2*1], [3/1]]", "array(array(+(1,1),*(2,1)),array(/(3,1)))")); // should fail later in conversion of AST to Arrays
	//EXPECT(match("[1, 2, 3] + sin([xyz])", ""));
	EXPECT(fail("[(]]", ""));
	//EXPECT(fail("(1)[1]", ""));
	//EXPECT(fail("xy[1]", ""));
	//EXPECT(fail("[1,2][1]", ""));
	//EXPECT(fail("[1,x[0]]", ""));
	//EXPECT(match("[1,2][0]", ""));

}


void test_subscription() {
	std::cout << "\ntest_subscription" << "\n";
	// simple indices
	EXPECT(match("a[1]", "[](`a`,1)"));
	EXPECT(fail("a[1][2]", "Parsing failed at: [2]"));
	EXPECT(match("a[1, 2]", "[](`a`,1,2)"));
	EXPECT(match("a[1, 2, 3]", "[](`a`,1,2,3)"));
	EXPECT(match("a[None]", "[](`a`,None(0))"));
	EXPECT(match("a[None, -1]", "[](`a`,None(0),-1)"));

	// slices
	EXPECT(match("a[0:1:-1]", "[](`a`,slice(0,1,-1))"));
	EXPECT(match("a[:]", "[](`a`,slice(None(0),None(0),None(0)))"));
	EXPECT(match("a[:1]", "[](`a`,slice(None(0),1,None(0)))"));
	EXPECT(match("a[1:]", "[](`a`,slice(1,None(0),None(0)))"));
    EXPECT(match("a[:-1:1]", "[](`a`,slice(None(0),-1,1))"));
	// EXPECT(match("a[::]", "[](`a`_,(None(0)_slice(None(0)_None(0)_None(0)))))")); // forbidden due to problems with grammar
	EXPECT(fail("a[::]", "Expected \"]\" at \":]\""));
	EXPECT(fail("a[:::]", "Expected \"]\" at \"::]\""));

	// multiple slices
	EXPECT(match("a[:, 1]", "[](`a`,slice(None(0),None(0),None(0)),1)"));
	EXPECT(match("a[:, None]", "[](`a`,slice(None(0),None(0),None(0)),None(0))"));

	//array index
	EXPECT(match("a[[1,3]]", "[](`a`,idxarray(1,3))"));
	EXPECT(match("a[[1,3], :3:-1]", "[](`a`,idxarray(1,3),slice(None(0),3,-1))"));
}


void test_operators() {
	std::cout << "\ntest_operators" << "\n";
	EXPECT(match("-1.23e-2", "-(0.0123)" ));
	EXPECT(match("2 ** 5", "**(2,5)"));
	EXPECT(match("2 ** 3 ** 4", "**(2,**(3,4))"));
	EXPECT(match("2 * 5", "*(2,5)"));
	EXPECT(match("-2 * 3 / 4", "/(*(-(2),3),4)"));
	EXPECT(match("1 * 5 / 4 + 3 - 4", "-(+(/(*(1,5),4),3),4)"));
	EXPECT(match("4 + 2 * (3 - 4)", "+(4,*(2,-(3,4)))"));
	EXPECT(match("((123))", "123"));
	EXPECT(match("((123)*1)/4", "/(*(123,1),4)"));
	EXPECT(match("(1*2) / (3*4)", "/(*(1,2),*(3,4))"));

	EXPECT(match("2 * -5", "*(2,-(5))")); // should possibly fail
	//ASSERT(! match("-1.23e-2 * -5.3"));
	//ASSERT(! match("+ -1.23e-2"));
	EXPECT(match("+ + -1.23e-2", "+(+(-(0.0123)))"));

	// relational operators
	EXPECT(match("2 > 1", ">(2,1)"));
	EXPECT(match("2 == 1", "==(2,1)"));
	EXPECT(match("2 <= 1", "<=(2,1)"));

	// boolean
	EXPECT(match("not 2 <= 1", "not(<=(2,1))"));
	EXPECT(match("not 2 <= 1 or False", "or(not(<=(2,1)),False(0))"));
	EXPECT(match("not 2 <= 1 or False and 1", "or(not(<=(2,1)),and(False(0),1))"));
	EXPECT(match("10 if True else 20", "ifelse(10,True(0),20)"));
}



int main()
{
	test_primary();
	test_arrays();
	test_subscription();
	test_operators();

}


#endif /* TEST_TEST_PARSER_CC_ */
