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
#include <memory>
#include "config.hh"
#include "assert.hh"
#include "VCL_v2_include.hh"
#include "arena_alloc.hh"


namespace bparser {
namespace details {


enum ResultStorage {
	none = 0,
	constant = 1,
	value = 2,
	temporary = 3,
	expr_result = 4,
	value_copy = 5
};

struct ScalarNode;
typedef std::shared_ptr<ScalarNode> ScalarNodePtr;

const int64_t true_value = 0xFFFFFFFFFFFFFFFFLL;
const int64_t false_value = 0x0000000000000000LL;

template<typename VecType>
void printVector(const VecType & v, const char * prefix)
{
    bool first = true;
    std::cout << prefix << "(";
    for(int i = 0; i < VecType::size(); i++)
    {
        if (first)
        {
            std::cout << v[i];
            first = false;
            continue;
        }

        std::cout << " ; " << v[i];
    }
    std::cout << ")" << std::endl;
}

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
	ScalarNodePtr  inputs_[3];
	// Number of (yet unprocessed) nodes depending on stored result.
	// Used in Processor to reuse temporary result storage.
	uint n_dep_nodes_;
	// index of the result in workspace, can be reused
	int result_idx_;

	unsigned char op_code_;
	std::string op_name_;
	// Pointer to (user provided) vector of values.
	double * values_;

	/**
	 * Factory functions fro special nodes.
	 */
	inline static ScalarNodePtr create_zero();
	inline static ScalarNodePtr create_one();
	inline static ScalarNodePtr create_const(double a);
	inline static ScalarNodePtr create_value(double *a);
	inline static ScalarNodePtr create_val_copy(double *a);
	inline static ScalarNodePtr create_result(ScalarNodePtr result, double *a);
	inline static ScalarNodePtr create_ifelse(ScalarNodePtr a, ScalarNodePtr b, ScalarNodePtr c);

	/**
	 * Generic factory functions for operation nodes.
	 * Separate template for different number of inputs (incoming edges in the graph).
	 */
	template <class T>
	static ScalarNodePtr  create(ScalarNodePtr a);

	template <class T>
	static ScalarNodePtr  create(ScalarNodePtr a, ScalarNodePtr b);



	ScalarNode()
	: result_storage(temporary),
	  n_inputs_(0),
	  n_dep_nodes_(0),
	  result_idx_(-1),
	  op_code_((unsigned char)0xff),
	  op_name_("none"),
	  values_(nullptr)
	{}

	virtual ~ScalarNode() {
	}

	void add_input(ScalarNodePtr  in)
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

	~ConstantNode() override {
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

	~ValueNode() override {
	}

};


struct ValueCopyNode : public ScalarNode {
	ValueCopyNode(double *ptr)
	{
		op_name_ = "ValueCopy";
		source_ptr_ = ptr;
		result_storage = value_copy;
		values_ = nullptr; // Is set automatically by Processor.
	}

	~ValueCopyNode() override {
	}

	double * source_ptr_;               ///< Pointer to data passed in constructor
};


struct ResultNode : public ScalarNode {
	ResultNode(ScalarNodePtr result_node, double *storage)
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

template<typename bool_type> struct b_to_d;

template<>
struct b_to_d<int64_t> {
    typedef double double_type;
};

template<>
struct b_to_d<Vec2db> {
    typedef Vec2d double_type;
};

template<>
struct b_to_d<Vec4db> {
    typedef Vec4d double_type;
};

template<>
struct b_to_d<Vec8db> {
    typedef Vec8d double_type;
};

template<typename bool_type> union b_to_d_mask;

template<>
union b_to_d_mask<int64_t> {
	int64_t	mask;
	double  value;
};

template<>
union b_to_d_mask<Vec2db> {
	Vec2db	mask;
	Vec2d  value;
};

template<>
union b_to_d_mask<Vec4db> {
	Vec4db	mask;
	Vec4d  value;
};

template<>
union b_to_d_mask<Vec8db> {
	Vec8db	mask;
	Vec8d  value;
};

template<typename bool_type>
inline typename b_to_d<bool_type>::double_type as_double(bool_type in) {
    b_to_d_mask<bool_type> x = {in};
    return x.value;
}

template<typename double_type> struct d_to_b;

template<>
struct d_to_b<double> {
    typedef int64_t bool_type;
};

template<>
struct d_to_b<Vec2d> {
    typedef Vec2db bool_type;
};

template<>
struct d_to_b<Vec4d> {
    typedef Vec4db bool_type;
};

template<>
struct d_to_b<Vec8d> {
    typedef Vec8db bool_type;
};

template<typename double_type> union d_to_b_mask;

template<>
union d_to_b_mask<double> {
	double  value;
	int64_t	mask;
};

template<>
union d_to_b_mask<Vec2d> {
	Vec2d  value;
	Vec2db	mask;
};

template<>
union d_to_b_mask<Vec4d> {
	Vec4d  value;
	Vec4db	mask;
};

template<>
union d_to_b_mask<Vec8d> {
	Vec8d  value;
	Vec8db	mask;
};


template<typename double_type>
inline typename d_to_b<double_type>::bool_type as_bool(double_type in) {
    d_to_b_mask<double_type> x = {in};
    return x.mask;
}

// previous functions

// union MaskDouble {
// 	int64_t	mask;
// 	double  value;
// };

// union DoubleMask {
// 	double  value;
// 	int64_t	mask;
// };

// inline double mask_to_double(int64_t x) {
// 	return reinterpret_cast<double &>(x);
// }

// inline int64_t double_to_mask(double x) {
// 	return reinterpret_cast<int64_t &>(x);
// }

// inline constexpr double mask_to_double(int64_t x) {
// 	MaskDouble m = {x};
// 	return m.value;
// }

// inline constexpr int64_t double_to_mask(double x) {
// 	DoubleMask m = {x};
// 	return m.mask;
// }

inline int64_t bitmask_false() {
	return false_value;
}

inline int64_t bitmask_true() {
	return true_value;
}

inline double double_false() {
	return as_double(bitmask_false());
}

inline double double_true() {
	return as_double(bitmask_true());
}

inline double double_bool(bool x) {
	return x ? double_true() : double_false();
}


/***
 * Declaration of particular scalar nodes for distinct operations.
 * Unique 'op_code' must be set to every particular operation ScalarNode.
 */

#define UNARY_FN(NAME, OP_CODE, FN) 									\
	struct NAME : public ScalarNode {									\
		static const char op_code = OP_CODE;							\
		static const char n_eval_args = 2;								\
		template <typename VecType>										\
		inline static void eval(VecType &res, VecType a) {				\
			res = FN(a);												\
		}																\
	}

struct _minus_ : public ScalarNode {
	static const char op_code = 1;
	static const char n_eval_args = 2;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a) {
		res = -a;
	}
};

struct _add_ : public ScalarNode {
	static const char op_code = 2;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		// std::cout << a << " + " << b << "\n";
		res = a + b;
	}
};

struct _sub_ : public ScalarNode {
	static const char op_code = 3;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		// std::cout << a << " - " << b << "\n";
 		res = a - b;
	}
};


struct _mul_ : public ScalarNode {
	static const char op_code = 4;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		// std::cout << a << " * " << b << "\n";
		res = a * b;
	}
};

struct _div_ : public ScalarNode {
	static const char op_code = 5;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		res = a / b;
	}
};

struct _mod_ : public ScalarNode {
	static const char op_code = 6;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _mod_::eval(VecType &res, VecType a, VecType b) {
	res = a - b * truncate(a / b);
}
template<>
inline void _mod_::eval<double>(double &res, double a, double b) {
	res = std::fmod(a, b);
}

struct _eq_ : public ScalarNode {
	static const char op_code = 7;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _eq_::eval(VecType &res, VecType a, VecType b) {
	res = as_double(a == b);
}
template<>
inline void _eq_::eval<double>(double &res, double a, double b) {
	res = double_bool(a == b);
}

struct _ne_ : public ScalarNode {
	static const char op_code = 8;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _ne_::eval(VecType &res, VecType a, VecType b) {
	res = as_double(a != b);
}
template<>
inline void _ne_::eval<double>(double &res, double a, double b) {
	res = double_bool(a != b);
}

struct _lt_ : public ScalarNode {
	static const char op_code = 9;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _lt_::eval(VecType &res, VecType a, VecType b) {
	std::cout << "In lt: " << std::endl;
	printVector<VecType>(a, "a");
	printVector<VecType>(b, "b");

	res = as_double(a < b);

	std::cout << "res pointer: " << &res << std::endl;

	printVector<VecType>(res, "res");
}
template<>
inline void _lt_::eval<double>(double &res, double a, double b) {
	res = double_bool(a < b);
}

struct _le_ : public ScalarNode {
	static const char op_code = 10;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _le_::eval(VecType &res, VecType a, VecType b) {
	res = as_double(a <= b);
}
template<>
inline void _le_::eval<double>(double &res, double a, double b) {
	res = double_bool(a <= b);
}

struct _neg_ : public ScalarNode {
	static const char op_code = 11;
	static const char n_eval_args = 2;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a);
};
template<typename VecType>
inline void _neg_::eval(VecType &res, VecType a) {
	res = as_double(!a);
}
template<>
inline void _neg_::eval<double>(double &res, double a) {
	res =  (a == double_true()) ? double_false() : double_true();	// we use bit masks for bool values
}

struct _or_ : public ScalarNode {
	static const char op_code = 12;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		res = as_double(as_bool(a) | as_bool(b));	// we use bit masks for bool values
	}
};

struct _and_ : public ScalarNode {
	static const char op_code = 13;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		res = as_double(as_bool(a) & as_bool(b));	// we use bit masks for bool values
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
	template <typename VecType>
	inline static void eval(VecType &res, VecType a);
};
template<typename VecType>
inline void _isnan_::eval(VecType &res, VecType a) {
	res = as_double(is_nan(a));
}
template<>
inline void _isnan_::eval<double>(double &res, double a) {
	res = double_bool(std::isnan(a));
}

struct _isinf_ : public ScalarNode {
	static const char op_code = 37;
	static const char n_eval_args = 2;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a);
};
template<typename VecType>
inline void _isinf_::eval(VecType &res, VecType a) {
	res = as_double(is_inf(a));
}
template<>
inline void _isinf_::eval<double>(double &res, double a) {
	res = double_bool(std::isinf(a));
}

struct _sgn_ : public ScalarNode {
	static const char op_code = 38;
	static const char n_eval_args = 2;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a);
};
template<typename VecType>
inline void _sgn_::eval(VecType &res, VecType a) {
	res = as_double(sign_bit(a));
}
template<>
inline void _sgn_::eval<double>(double &res, double a) {
	res =  a > 0 ? 1.0 : (a < 0 ? -1.0 : 0.0);
}

struct _atan2_ : public ScalarNode {
	static const char op_code = 39;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		// TODO: vectorize
		res = atan2(a, b);
	}
};

struct _pow_ : public ScalarNode {
	static const char op_code = 40;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b) {
		// TODO: vectorize
		res = pow(a, b);
	}
};

struct _max_ : public ScalarNode {
	static const char op_code = 41;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _max_::eval(VecType &res, VecType a, VecType b) {
	res = max(a, b);
}
template<>
inline void _max_::eval<double>(double &res, double a, double b) {
	//std::cout << "max " << a << "," << b << "\n";
	res = (a > b) ? a : b;
}

struct _min_ : public ScalarNode {
	static const char op_code = 42;
	static const char n_eval_args = 3;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b);
};
template<typename VecType>
inline void _min_::eval(VecType &res, VecType a, VecType b) {
	res = min(a, b);
}
template<>
inline void _min_::eval<double>(double &res, double a, double b) {
	//std::cout << "min " << a << "," << b << "\n";
	res = (a > b) ? b : a;
}

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
	template <typename VecType>
	inline static void eval(VecType &res, VecType a) {
		res = a;
		//std::cout << a << " -copy-> " << res << "\n";
	}
};


struct _ifelse_ : public ScalarNode {
	static const char op_code = 51;
	static const char n_eval_args = 4;
	template <typename VecType>
	inline static void eval(VecType &res, VecType a, VecType b, VecType c);
};
template<typename VecType>
inline void _ifelse_::eval(VecType &res, VecType a, VecType b, VecType c) {
	std::cout << "In select: " << std::endl;
	printVector<VecType>(a, "a");
	printVector<VecType>(b, "b");
	printVector<VecType>(c, "c");
	
	res = select(as_bool(b), a, c);	// we use bit masks for bool values

	// printVector<d_to_b<VecType>::bool_type>(as_bool(b), "b_bool");
	printVector<VecType>(res, "res");
}
template<>
inline void _ifelse_::eval<double>(double &res, double a, double b, double c) {
	res = as_bool(b) ? a : c;		// we use bit masks for bool values
}
UNARY_FN(_log2_, 	52, log2);

/***********************
 * Construction Nodes.
 */

ScalarNodePtr ScalarNode::create_zero() {
	static ScalarNodePtr zero = std::make_shared<ConstantNode>(0.0);
	return zero;
}

ScalarNodePtr ScalarNode::create_one() {
	static ScalarNodePtr one = std::make_shared<ConstantNode>(1.0);
	return one;
}


inline ScalarNodePtr ScalarNode::create_const(double a) {
	return std::make_shared<ConstantNode>(a);
}

// create value node
inline ScalarNodePtr ScalarNode::create_value(double *a)  {
	return std::make_shared<ValueNode>(a);
}

// create value node
inline ScalarNodePtr ScalarNode::create_val_copy(double *a)  {
	return std::make_shared<ValueCopyNode>(a);
}

// create result node
inline ScalarNodePtr ScalarNode::create_result(ScalarNodePtr result, double *a)  {
	BP_ASSERT(result->result_storage != none);
	if (result->result_storage != temporary) {
		result = ScalarNode::create<_copy_>(result);
	}
	result->values_ = a;
	result->result_storage = expr_result;
	return result;
}

template <class T>
ScalarNodePtr ScalarNode::create(ScalarNodePtr a) {
	std::shared_ptr<T> node_ptr = std::make_shared<T>();
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

	return node_ptr;
}

template <class T>
ScalarNodePtr ScalarNode::create(ScalarNodePtr a, ScalarNodePtr b) {
	std::shared_ptr<T> node_ptr = std::make_shared<T>();
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

	return node_ptr;
}

inline ScalarNodePtr  ScalarNode::create_ifelse(ScalarNodePtr a, ScalarNodePtr b, ScalarNodePtr c)  {
	std::shared_ptr<_ifelse_> node_ptr = std::make_shared<_ifelse_>();
	node_ptr->op_code_ = _ifelse_::op_code;
	node_ptr->set_name(typeid(_ifelse_).name());
	node_ptr->add_input(a);
	node_ptr->add_input(b);
	node_ptr->add_input(c);
	BP_ASSERT(_ifelse_::n_eval_args == 4);
	node_ptr->result_storage = temporary;

	return node_ptr;
}



} // namespace details
} // namespace bparser
#endif /* INCLUDE_SCALAR_NODE_HH_ */
