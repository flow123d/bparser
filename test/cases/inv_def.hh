#ifndef NITPICK_IDE_IGNORE
#include "nitpick_include.hh"
#endif //NITPICK_IDE_IGNORE

void def(ExprCase& c) {

	// parse an expression.
	c.parse("inv(m3)");

	c.set_variable("m3", { 3,3 });

	// Set the result variable shape
	c.set_result_shape({3,3});


}