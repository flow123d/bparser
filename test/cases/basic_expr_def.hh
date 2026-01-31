/*
* Prepare an environment for both DAG generation and then subsequent running
*/

#ifndef NITPICK_IDE_IGNORE
#include "nitpick_include.hh"
#endif //NITPICK_IDE_IGNORE

void def(ExprCase& c) {

	// parse an expression.
	c.parse("1 * v1 + cs1 * v2");

	// "cs1" constant with shape {}, i.e. scalar and values {2}.
	c.set_constant("cs1", {}, {2});
	// "cv1" vector constant with shape {3}
	//P_SET_CONSTANT(cv1, {3}, ARG({1, 2, 3}));
	// "v1" variable with shape {3}; v1 is pointer to the value space
	c.set_variable("v1", {3});
	c.set_variable("v2", {3});
	// Set the result variable (the last line of the expression)
	c.set_result_shape({ 3 });


}
