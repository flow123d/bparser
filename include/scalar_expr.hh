/*
 * scalar_expr.hh
 *
 *  Created on: Jan 4, 2020
 *      Author: jb
 */

#ifndef INCLUDE_SCALAR_EXPR_HH_
#define INCLUDE_SCALAR_EXPR_HH_

#include <vector>
#include <cmath>
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
	char op_code_;
	double * values_;

	/**
	 * Factory functions fro special nodes.
	 */
	static ScalarNode * create_const(double a);
	static ScalarNode * create_value(double *a);
	static ScalarNode * create_result(ScalarNode *result, double *a);

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
	  op_code_(0xff),
	  values_(nullptr)
	{}

	void add_input(ScalarNode * in)
	{
		ASSERT(n_inputs_ < 3);
		inputs_[n_inputs_] = in;
		n_inputs_+=1;
	}

	double * get_value() {
		ASSERT(values_ != nullptr);
		return values_;
	}
};




struct ConstantNode : public ScalarNode {
	ConstantNode(double v)
	: value_(v)
	{
		values_ = &value_;
		result_storage = constant;
	}

	double value_;
};


struct ValueNode : public ScalarNode {
	ValueNode(double *ptr)
	{
		values_ = ptr;
		result_storage = value;
	}
};


struct ResultNode : public ScalarNode {
	ResultNode(ScalarNode *result_node, double *storage)
	{
		values_ = storage;
		ScalarNode::add_input(result_node);
		result_storage = expr_result; // use slot of the result of value, i.e. inputs[0]
	}
};




/***********************
 * Operation Nodes.
 */

double mask_to_double(int64_t x) {
	return *(reinterpret_cast<const double *>(& x));
}

int64_t double_to_mask(double x) {
	return *(reinterpret_cast<const int64_t *>(& x));
}

static const int64_t bitmask_false = 0x0000000000000000L;
static const int64_t bitmask_true = 0x1111111111111111L;
static const double double_false = mask_to_double(bitmask_false);
static const double double_true = mask_to_double(bitmask_true);


double double_bool(bool x) {
	return x ? double_true : double_false;
}

#define UNARY_FN(NAME, OP_CODE, FN) 									\
	struct NAME : public ScalarNode {									\
		static const char op_code = OP_CODE;							\
		static const char n_eval_args = 2;								\
		inline static void eval(double &res, double a) {				\
			res = FN(a);												\
		}																\
	};

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
		res = a + b;
	}
};

struct _sub_ : public ScalarNode {
	static const char op_code = 3;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		//std::cout << a << " - " << b << "\n";
 		res = a - b;
	}
};


struct _mul_ : public ScalarNode {
	static const char op_code = 4;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
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
		res =  a == double_true ? double_false : double_true;		// we use bit masks for bool values
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
	static const char op_code = 40;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res = a;
		//std::cout << a << " -copy-> " << res << "\n";
	}
};




/***********************
 * Construction Nodes.
 */



ScalarNode * ScalarNode::create_const(double a) {
	return new ConstantNode(a);
//		nodes.push_back(node);
//		n_constants += 1;
//		return node;
}

// create value node
ScalarNode * ScalarNode::create_value(double *a)  {
	return new ValueNode(a);
//		nodes.push_back(node);
//		n_values += 1;
//		return node;
}

// create result node
ScalarNode * ScalarNode::create_result(ScalarNode *result, double *a)  {
	if (result->result_storage == constant || result->result_storage == value) {
		result = ScalarNode::create<_copy_>(result);
	}

	return new ResultNode(result, a);
//		nodes.push_back(node);
//		results.push_back(node);
//		return node;
}

template <class T>
ScalarNode * ScalarNode::create(ScalarNode *a) {
	T * node_ptr = new T();
	node_ptr->op_code_ = T::op_code;
	node_ptr->add_input(a);
	if (T::n_eval_args == 1) {
		node_ptr->result_storage = none;
	} else {
		ASSERT(T::n_eval_args == 2);
		node_ptr->result_storage = temporary;
	}

//		nodes.push_back(node_ptr);
	return node_ptr;
}

template <class T>
ScalarNode * ScalarNode::create(ScalarNode *a, ScalarNode *b) {
	T * node_ptr = new T();
	node_ptr->op_code_ = T::op_code;
	node_ptr->add_input(a);
	node_ptr->add_input(b);
	if (T::n_eval_args == 2) {
		node_ptr->result_storage = none;
	} else {
		ASSERT(T::n_eval_args == 3);
		node_ptr->result_storage = temporary;
	}

//		nodes.push_back(node_ptr);
	return node_ptr;
}


/**
 * Auxiliary class to form evaluation graph and make topological sort
 * to get order of operations.
 * TODO: optimize topological sort for the number of temporaries.
 */
struct ScalarExpression {
	typedef std::vector<ScalarNode *> NodeVec;

	NodeVec nodes;
	NodeVec sorted;
	NodeVec results;
	uint n_constants;
	uint n_values;
	std::vector<uint> storage;

	ScalarExpression(std::vector<ScalarNode *> res)
	:
		results(res.begin(), res.end()),
		n_constants(0),
		n_values(0)
	{
		// ScalarNode::reslut_idx_ == -1,
		// we set it to -2 to identify passed nodes
		nodes.clear();
		for(ScalarNode * node : results)
			add_node(node);

		for(uint i=0; i < nodes.size(); ++i) {
			ScalarNode * node = nodes[i];
			for(uint in=0; in < node->n_inputs_; ++in)  {
				ScalarNode * other = node->inputs_[in];
				if (other->result_idx_ != -2) {
					ASSERT(other->result_idx_ == -1);
					add_node(other);
				}
			}


		}

		ASSERT(sorted.size() == 0);
		sort_nodes();
	}

	void add_node(ScalarNode * node) {
		node->result_idx_ = -2;
		nodes.push_back(node);
		if (node->result_storage == constant)
			n_constants++;
		if (node->result_storage == value)
			n_values++;

	}


	~ScalarExpression() {
		for(ScalarNode *node : nodes) {
			delete node;
		}
	}


	// Perform topological sort, set temporary indices
	NodeVec & sort_nodes() {
        if (sorted.size() > 0)
        	return sorted;

		for(ScalarNode *node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNode *node : nodes)
			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ += 1;
				//std::cout << "  node: " << node->inputs_[in] << " n_dep: " << node->inputs_[in]->n_dep_nodes_ << "\n";
			}

		NodeVec stack(results.begin(), results.end());

		while (stack.size() > 0) {
			ScalarNode * node = stack.back();
			//std::cout << "node: " << node << " res: " << node->result_storage << "\n";
			stack.pop_back();
			sorted.push_back(node);

			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ -= 1;
				//std::cout << "  node: " << node->inputs_[in] << " n_dep: " << node->inputs_[in]->n_dep_nodes_ << "\n";
				if (node->inputs_[in]->n_dep_nodes_ == 0)
					stack.push_back(node->inputs_[in]);
			}
		}

		setup_result_storage();
		return sorted;
	}

	void setup_result_storage() {
		for(ScalarNode *node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNode *node : nodes)
			for(uint in=0; in < node->n_inputs_; ++in) node->inputs_[in]->n_dep_nodes_ += 1;

		for(auto it=sorted.rbegin(); it != sorted.rend(); ++it) {
			ScalarNode * node = *it;
			allocate_storage(node);
			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ -= 1;
				if (node->inputs_[in]->n_dep_nodes_ == 0) {
					deallocate_storage(node->inputs_[in]);
				}
			}
		}
	}
	/**
	 * forward processing in topological order
	 * for all nodes set N.n_dep
	 * allocate temporary for node N if: all N->input[i]  are processed
	 * if ++(M=N->input[i]).n_dep == M.max_dep deallocate temporary of M
	 *
	 * <=>
	 *
	 * backward processing ??
	 * allocate if all inputs are unprocessed ... first triggered
	 * deallocate M if --M.n_dep == 0
	 *
	 */
	void allocate_storage(ScalarNode *node) {

		if (node->result_storage == none || node->result_storage == expr_result)
		{
			ScalarNode * reused_node = node->inputs_[0];
			if (reused_node->result_storage == temporary || reused_node->result_storage == none) {
				// reuse previous slot
				node->result_idx_ = reused_node->result_idx_;
			} else {
				ASSERT(false);
				// result and operations with none storage
				// can not follow ConstantNode or ValueNode
			}
		} else {
			//std::cout << "node: " << node << " get slot: " << storage.size() << "\n";
			node->result_idx_ = get_free_slot();
		}
	}

	uint get_free_slot() {
		for(uint i=0; i<storage.size(); ++i)
			if (storage[i] == 0) {
				storage[i] = 1;
				return i;
			}
		storage.push_back(1);
		return storage.size() - 1;
	}

	void deallocate_storage(ScalarNode *node) {
		if (node->result_storage == temporary) {
			//std::cout << "dealoc node: " << node << " idx: " << node->result_idx_ << "\n";
			storage[node->result_idx_] = 0;
		}
	}


	uint n_vectors() {
		return storage.size();
	}

};






} // namespace details
} // namespace bparser
#endif /* INCLUDE_SCALAR_EXPR_HH_ */
