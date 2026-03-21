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
#include "nitpick_common.hh"

//void def(ExprCase c);
//ExpressionDAG gen(const NodeMap& node_map);



enum VariableAllocationEnum {
	KeepDefault,
	ForceVariable,
	ForceVarCopy,
	ForceConst
};

VariableAllocationEnum globalVariableAllocationMethod = KeepDefault;

ScalarNodePtr create_variable_node(double* ptr, VarType preffered_type) {
	using namespace bparser::details;
	switch (globalVariableAllocationMethod) {
	case KeepDefault: {
		switch (preffered_type) {
		case Variable:
			return ScalarNode::create_value(ptr);
		case VarCopy:
			return ScalarNode::create_val_copy(ptr);
		case Const:
			return ScalarNode::create_const(*ptr);
		case ConstBool:
			return ScalarNode::create_const_bool(*ptr);
		}
	}
	case ForceVariable:
		return ScalarNode::create_value(ptr);
	case ForceVarCopy:
		return ScalarNode::create_val_copy(ptr);
	case ForceConst:
		return ScalarNode::create_const(*ptr);
	}
}


#endif //NITPICK_INCLUDE_HH