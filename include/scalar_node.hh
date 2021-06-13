/*
 * scalar_expr.hh
 *
 *  Created on: Jan 4, 2020
 *      Author: jb
 */

#ifndef INCLUDE_SCALAR_NODE_HH_
#define INCLUDE_SCALAR_NODE_HH_

#include <vector>
#include <cmath>
#include <map>
#include <typeinfo>
#include "config.hh"
#include "assert.hh"



namespace bparser {
namespace details {


enum ResultStorage {
	none = 0,
	constant = 1,
	value = 2,
	temporary = 3,
	expr_result = 4
};

/**
 * ScalarNodes describes DAG of elementary operations. From this
 * description of the expression the Processor instance is constructed
 * for efficient evaluation of the expression.
 *
 * ScalarNodes are not meant to be directly used.
 * Use Array class to construct general vector expressions.
 */

struct ScalarNode {
	static const char terminate_op_code = 0;

	ResultStorage result_storage;
	uint n_inputs_;
	ScalarNode * inputs_[3];
	// Number of (yet unprocessed) nodes depending on stored result.
	// Used in Processor to reuse temporary result storage.
	uint n_dep_nodes_;
	// index of the result in workspace, can be reused
	int result_idx_;
	unsigned char op_code_;
	std::string op_name_;
	double * values_;

	/**
	 * Factory functions fro special nodes.
	 */
	static ScalarNode * create_const(double a);
	static ScalarNode * create_value(double *a);
	static ScalarNode * create_result(ScalarNode *result, double *a);
	static ScalarNode * create_ifelse(ScalarNode *a, ScalarNode *b, ScalarNode *c);

	/**
	 * Generic factory functions for operation nodes.
	 * Separate template for different number of inputs (incoming edges in the graph).
	 */
	template <class T>
	static ScalarNode * create(ScalarNode *a);

	template <class T>
	static ScalarNode * create(ScalarNode *a, ScalarNode *b);



	ScalarNode()
	: result_storage(temporary),
	  n_inputs_(0),
	  n_dep_nodes_(0),
	  result_idx_(-1),
	  op_code_((unsigned char)0xff),
	  op_name_("none"),
	  values_(nullptr)
	{}

	void add_input(ScalarNode * in)
	{
		BP_ASSERT(n_inputs_ < 3);
		inputs_[n_inputs_] = in;
		n_inputs_+=1;
	}

	double * get_value() {
		BP_ASSERT(values_ != nullptr);
		return values_;
	}

	void set_name(std::string str) {
		uint a = str.find_first_of('_')+1;
		uint b = str.find_last_of('_');
		op_name_ = str.substr(a, b-a);
	}
};




struct ConstantNode : public ScalarNode {
	ConstantNode(double v)
	: value_(v)
	{
		op_name_ = "Const";
		values_ = &value_;
		result_storage = constant;
	}

	double value_;
};


struct ValueNode : public ScalarNode {
	ValueNode(double *ptr)
	{
		op_name_ = "Value";
		values_ = ptr;
		result_storage = value;
	}
};


struct ResultNode : public ScalarNode {
	ResultNode(ScalarNode *result_node, double *storage)
	{
		op_name_ = "Result";
		values_ = storage;
		ScalarNode::add_input(result_node);
		result_storage = expr_result; // use slot of the result of value, i.e. inputs[0]
	}
};




/***********************
 * Operation Nodes.
 */

union MaskDouble {
	int64_t	mask;
	double  value;
};

union DoubleMask {
	double  value;
	int64_t	mask;
};

inline double mask_to_double(int64_t x) {
    MaskDouble a = {x};
    return a.value;
    //return reinterpret_cast<double &>(x);
}

inline int64_t double_to_mask(double x) {
    DoubleMask a = {x};
    return a.mask;
}

//inline constexpr double mask_to_double(int64_t x) {
//	MaskDouble m = {x};
//	return m.value;
//}
//
//inline constexpr int64_t double_to_mask(double x) {
//	DoubleMask m = {x};
//	return m.mask;
//}

inline int64_t bitmask_false() {
	return 0x0000000000000000L;
}

inline int64_t bitmask_true() {
	return 0x1111111111111111L;
}

inline double double_false() {
	return mask_to_double(bitmask_false());
}

inline double double_true() {
	return mask_to_double(bitmask_true());
}


inline double double_bool(bool x) {
	return x ? double_true() : double_false();
}

#define UNARY_FN(NAME, OP_CODE, FN) 									\
	struct NAME : public ScalarNode {									\
		static const char op_code = OP_CODE;							\
		static const char n_eval_args = 2;								\
		inline static void eval(double &res, double a) {				\
			res = FN(a);												\
		}																\
	}

struct _minus_ : public ScalarNode {
	static const char op_code = 1;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res = -a;
	}
};


struct _add_ : public ScalarNode {
	static const char op_code = 2;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// std::cout << a << " + " << b << "\n";
		res = a + b;
	}
};

struct _sub_ : public ScalarNode {
	static const char op_code = 3;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// std::cout << a << " - " << b << "\n";
 		res = a - b;
	}
};


struct _mul_ : public ScalarNode {
	static const char op_code = 4;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// std::cout << a << " * " << b << "\n";
		res = a * b;
	}
};

struct _div_ : public ScalarNode {
	static const char op_code = 5;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		res = a / b;
	}
};

struct _mod_ : public ScalarNode {
	static const char op_code = 6;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = std::fmod(a, b);
	}
};

struct _eq_ : public ScalarNode {
	static const char op_code = 7;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = double_bool(a == b);
	}
};

struct _ne_ : public ScalarNode {
	static const char op_code = 8;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = double_bool(a != b);
	}
};


struct _lt_ : public ScalarNode {
	static const char op_code = 9;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = double_bool(a < b);
	}
};

struct _le_ : public ScalarNode {
	static const char op_code = 10;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = double_bool(a <= b);
	}
};

struct _neg_ : public ScalarNode {
	static const char op_code = 11;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		// TODO: vectorize
		res =  (a == double_true()) ? double_false() : double_true();		// we use bit masks for bool values
	}
};


struct _or_ : public ScalarNode {
	static const char op_code = 12;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = mask_to_double( double_to_mask(a) | double_to_mask(b));	// we use bit masks for bool values
	}
};

struct _and_ : public ScalarNode {
	static const char op_code = 13;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res = mask_to_double( double_to_mask(a) & double_to_mask(b));	// we use bit masks for bool values
	}
};


UNARY_FN(_abs_, 	20, abs);
UNARY_FN(_sqrt_, 	21, sqrt);
UNARY_FN(_exp_, 	22, exp);
UNARY_FN(_log_, 	23, log);
UNARY_FN(_log10_, 	24, log10);
UNARY_FN(_sin_, 	25, sin);
UNARY_FN(_sinh_, 	26, sinh);
UNARY_FN(_asin_, 	27, asin);
UNARY_FN(_cos_, 	28, cos);
UNARY_FN(_cosh_, 	29, cosh);
UNARY_FN(_acos_, 	30, acos);
UNARY_FN(_tan_, 	31, tan);
UNARY_FN(_tanh_, 	32, tanh);
UNARY_FN(_atan_, 	33, atan);
UNARY_FN(_ceil_, 	34, ceil);
UNARY_FN(_floor_, 	35, floor);



struct _isnan_ : public ScalarNode {
	static const char op_code = 36;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		// TODO: vectorize
		res = double_bool(std::isnan(a));
	}
};

struct _isinf_ : public ScalarNode {
	static const char op_code = 37;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		// TODO: vectorize
		res = double_bool(std::isinf(a));
	}
};

struct _sgn_ : public ScalarNode {
	static const char op_code = 38;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		// TODO: vectorize
		res =  a > 0 ? 1.0 : (a < 0 ? -1.0 : 0.0);
	}
};

struct _atan2_ : public ScalarNode {
	static const char op_code = 39;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res =  atan2(a, b);
	}
};

struct _pow_ : public ScalarNode {
	static const char op_code = 40;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		res =  pow(a, b);
	}
};

struct _max_ : public ScalarNode {
	static const char op_code = 41;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		//std::cout << "max " << a << "," << b << "\n";
		res =  (a>b) ? a : b;
	}
};

struct _min_ : public ScalarNode {
	static const char op_code = 42;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		// TODO: vectorize
		//std::cout << "min " << a << "," << b << "\n";
		res =  (a>b) ? b : a;
	}
};

//struct _iadd_ : public ScalarNode {
//	static const char op_code = 2;
//	static const char n_eval_args = 2;
//	inline static void eval(double &res, double a) {
//		res += a;
//		//std::cout << a << " -iadd-> " << res << "\n";
//	}
//};
//
//struct _imul_ : public ScalarNode {
//	static const char op_code = 4;
//	static const char n_eval_args = 2;
//	inline static void eval(double &res, double a) {
//		res *= a;
//		//std::cout << a << " -iadd-> " << res << "\n";
//	}
//};


// Explicit copy is used only for expression result and only in the case that there is no
// operation between inputs and the result.
struct _copy_ : public ScalarNode {
	static const char op_code = 50;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res = a;
		//std::cout << a << " -copy-> " << res << "\n";
	}
};


struct _ifelse_ : public ScalarNode {
	static const char op_code = 51;
	static const char n_eval_args = 4;
	inline static void eval(double &res, double a, double b, double c) {
		// TODO: vectorize
		res = double_to_mask(b) ? a : c;	// we use bit masks for bool values
	}
};




/***********************
 * Construction Nodes.
 */



inline ScalarNode * ScalarNode::create_const(double a) {
	return new ConstantNode(a);
//		nodes.push_back(node);
//		n_constants += 1;
//		return node;
}

// create value node
inline ScalarNode * ScalarNode::create_value(double *a)  {
	return new ValueNode(a);
//		nodes.push_back(node);
//		n_values += 1;
//		return node;
}

// create result node
inline ScalarNode * ScalarNode::create_result(ScalarNode *result, double *a)  {
	BP_ASSERT(result->result_storage != none);
	if (result->result_storage != temporary) {
		result = ScalarNode::create<_copy_>(result);
	}
	result->values_ = a;
	result->result_storage = expr_result;
	return result;
}

template <class T>
ScalarNode * ScalarNode::create(ScalarNode *a) {
	T * node_ptr = new T();
	node_ptr->op_code_ = T::op_code;
	node_ptr->set_name(typeid(T).name());
	node_ptr->add_input(a);
	if (T::n_eval_args == 1) {
		// Note: in place operations are not supported
		node_ptr->result_storage = none;
	} else {
		BP_ASSERT(T::n_eval_args == 2);
		node_ptr->result_storage = temporary;
	}

//		nodes.push_back(node_ptr);
	return node_ptr;
}

template <class T>
ScalarNode * ScalarNode::create(ScalarNode *a, ScalarNode *b) {
	T * node_ptr = new T();
	node_ptr->op_code_ = T::op_code;
	node_ptr->set_name(typeid(T).name());
	node_ptr->add_input(a);
	node_ptr->add_input(b);
	if (T::n_eval_args == 2) {
		// Note: in place operations are not supported
		node_ptr->result_storage = none;
	} else {
		BP_ASSERT(T::n_eval_args == 3);
		node_ptr->result_storage = temporary;
	}

//		nodes.push_back(node_ptr);
	return node_ptr;
}

inline ScalarNode * ScalarNode::create_ifelse(ScalarNode *a, ScalarNode *b, ScalarNode *c)  {
	auto * node_ptr = new _ifelse_;
	node_ptr->op_code_ = _ifelse_::op_code;
	node_ptr->set_name(typeid(_ifelse_).name());
	node_ptr->add_input(a);
	node_ptr->add_input(b);
	node_ptr->add_input(c);
	BP_ASSERT(_ifelse_::n_eval_args == 4);
	node_ptr->result_storage = temporary;

//		nodes.push_back(node_ptr);
	return node_ptr;
}



} // namespace details
} // namespace bparser
#endif /* INCLUDE_SCALAR_NODE_HH_ */
