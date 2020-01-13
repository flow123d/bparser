/*
 * test_parser.cc
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */

#include <string>
#include "assert.hh"
#include "parser.hh"

void test_expr(std::string expr) {
	using namespace bparser;
	Parser p(4);
	p.parse(expr);
}


int main()
{
	test_expr("123.0");
}




