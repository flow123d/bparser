
#ifndef NITPICK_IDE_IGNORE
# include "nitpick/nitpick_include.hh"
#endif // NITPICK_IDE_IGNORE

void def(ExprCase& c) {

c.parse("norm2(v1)");

c.set_variable("v1", { 3 });
c.set_result_shape({ });
}

