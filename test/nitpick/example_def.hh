
#ifndef NITPICK_IDE_IGNORE
# include "nitpick_include.hh"
#endif //This block will not be compiled during generation/run. But it helps the IDE understand the code

void def(ExprCase& c) {

c.parse("1+1");

//c.set_variable("v1" {3});
c.set_result_shape({});
}