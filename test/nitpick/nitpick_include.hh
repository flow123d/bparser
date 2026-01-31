#ifndef NITPICK_INCLUDE_HH
#define NITPICK_INCLUDE_HH

#ifndef NITPICK_DEF_FILE
# ifdef DEF_FILE
#  define NITPICK_DEF_FILE DEF_FILE
# else
#  error [Nitpick] DEF_FILE not defined
# endif
#endif

#ifndef NITPICK_GEN_FILE
# ifdef GEN_FILE
#  define NITPICK_GEN_FILE GEN_FILE
# else
#  error [Nitpick] GEN_FILE not defined
# endif
#endif

#include "test_tools.hh"
#include "parser.hh"
#include <iostream>
#include "exprcase.hh"

//void def(ExprCase c);
//ExpressionDAG gen(const NodeMap& node_map);


#endif //NITPICK_INCLUDE_HH