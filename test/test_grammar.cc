/*
 * test_parser.hh
 *
 *  Created on: Dec 26, 2019
 *      Author: jb
 */

#ifndef TEST_TEST_PARSER_CC_
#define TEST_TEST_PARSER_CC_

#include <string>

#include "../include/grammar_old.hh"
#include "assert.hh"

using boost::spirit::ascii::space;
using boost::spirit::utree;

typedef std::string::const_iterator iterator_type;
typedef bparser::grammar<iterator_type> Grammar;


bool match(std::string s) {

    Grammar g; // Our grammar

    std::string s1 = s;
    std::string::const_iterator iter = s1.begin();
    std::string::const_iterator end = s1.end();

    bool result = phrase_parse(iter, end, g, space);
    //std::cout << "Expr: '" << s << "' Resutl: " << result << "\n";
    return result && iter == end;
}




int main()
{
	//bool r = match("-1.23e-2 * 5.3 / 4 + 3 - 4");
	bool r = match("+ + -1.23e-2");
	std::cout << r << "\n";


	// Test that grammar match valid expressions and fails for invalid.
	BP_ASSERT(match(("123")));
	BP_ASSERT(match(("123.0")));
	BP_ASSERT(match(("1.23e2")));
	BP_ASSERT(match(("1.23e-2")));
	BP_ASSERT(match(("-1.23e-2")));
	BP_ASSERT(match(("1.23e-2 * 5.3")));
	BP_ASSERT(match(("-1.23e-2 * 5.3 / 4")));
	BP_ASSERT(match(("-1.23e-2 * 5.3 / 4 + 3 - 4")));
	BP_ASSERT(match(("-1.23e-2 * 5.3 / 4 + 2 * (3 - 4)")));
	BP_ASSERT(match(("((123))")));
	BP_ASSERT(match(("((123)*1)/4")));
	BP_ASSERT(match(("(1*2) / (3*4)")));

	//ASSERT(! match(("-1.23e-2 * -5.3")));
	//ASSERT(! match(("-1.23e-2 * -5.3")));
	//ASSERT(! match(("+ -1.23e-2")));
	BP_ASSERT(! match(("+ + -1.23e-2")));
}


#endif /* TEST_TEST_PARSER_CC_ */
