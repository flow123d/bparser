#ifndef NITPICK_IDE_IGNORE
#include "nitpick_include.hh"
#endif //NITPICK_IDE_IGNORE

void def(ExprCase& c) {

	// parse an expression.
	c.parse("inv(m4)");

	c.set_variable("m4", { 4,4 });

	// Set the result variable shape
	c.set_result_shape({ 4,4 });


}