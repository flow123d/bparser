/*
 * test_processor.cc
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */



#include <string>
#include "assert.hh"
#include "processor.hh"




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
	ASSERT(match(("123")));
	ASSERT(match(("123.0")));
	ASSERT(match(("1.23e2")));
	ASSERT(match(("1.23e-2")));
	ASSERT(match(("-1.23e-2")));
	ASSERT(match(("1.23e-2 * 5.3")));
	ASSERT(match(("-1.23e-2 * 5.3 / 4")));
	ASSERT(match(("-1.23e-2 * 5.3 / 4 + 3 - 4")));
	ASSERT(match(("-1.23e-2 * 5.3 / 4 + 2 * (3 - 4)")));
	ASSERT(match(("((123))")));
	ASSERT(match(("((123)*1)/4")));
	ASSERT(match(("(1*2) / (3*4)")));

	//ASSERT(! match(("-1.23e-2 * -5.3")));
	//ASSERT(! match(("-1.23e-2 * -5.3")));
	//ASSERT(! match(("+ -1.23e-2")));
	ASSERT(! match(("+ + -1.23e-2")));
}

