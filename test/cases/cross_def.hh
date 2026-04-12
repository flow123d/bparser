
#ifndef NITPICK_IDE_IGNORE
# include "nitpick/nitpick_include.hh"
#endif // NITPICK_IDE_IGNORE

void def(ExprCase& c) {

	c.parse("cross(u,v)");

	c.set_variable("u", { 3 });
	c.set_variable("v", { 3 });
	c.set_result_shape( { 3 });
}

