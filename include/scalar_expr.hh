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


struct _add_ : public ScalarNode {
	static const char op_code = 1;
	static const char n_eval_args = 3;
	inline static void eval(double &res, double a, double b) {
		res = a + b;
	}
};

struct _iadd_ : public ScalarNode {
	static const char op_code = 2;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res += a;
	}
};


// Explicit copy is used only for expression result and only in the case that there is no
// operation between inputs and the result.
struct _copy_ : public ScalarNode {
	static const char op_code = 3;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res = a;
	}
};

struct _abs_ : public ScalarNode {
	static const char op_code = 4;
	static const char n_eval_args = 2;
	inline static void eval(double &res, double a) {
		res = fabs(a);
	}
};




struct ScalarExpression {
	typedef std::vector<ScalarNode *> NodeVec;

	NodeVec nodes;
	NodeVec sorted;
	NodeVec results;
	uint n_constants;
	uint n_values;
	uint next_result_idx;
	std::vector<uint> storage;

	ScalarExpression()
	: n_constants(0), n_values(0), next_result_idx(0)
	{}

	// create const node
	ScalarNode * create_const(double a) {
		auto node = new ConstantNode(a);
		nodes.push_back(node);
		n_constants += 1;
	}

	// create value node
	ScalarNode * create_value(double *a)  {
		auto node = new ValueNode(a);
		nodes.push_back(node);
		n_values += 1;
	}

	// create result node
	ScalarNode * create_result(ScalarNode *result, double *a)  {
		if (result->result_storage == constant || result->result_storage == value) {
			result = create<_copy_>(result);
		}

		auto node = new ResultNode(result, a);
		nodes.push_back(node);
		results.push_back(node);
	}

	template <class T>
	ScalarNode * create(ScalarNode *a);

	template <class T>
	ScalarNode * create(ScalarNode *a, ScalarNode *b);

	template <class T, uint n_eval_args>
	void check_eval();

	~ScalarExpression() {
		for(ScalarNode *node : nodes) {
			delete node;
		}
	}


	// Perform topological sort, set temporary indices
	NodeVec & sort_nodes() {
		for(ScalarNode *node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNode *node : nodes)
			for(uint in=0; in < node->n_inputs_; ++in) node->inputs_[in]->n_dep_nodes_ += 1;

		NodeVec stack(results.begin(), results.end());
		while (stack.size() > 0) {
			ScalarNode * node = stack.back();
			stack.pop_back();
			sorted.push_back(node);

			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ -= 1;
				if (node->inputs_[in]->n_dep_nodes_ == 0)
					stack.push_back(node);
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
					deallocate_storage(node);
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
			storage[node->result_idx_] = 0;
		}
	}


	uint n_vectors() {
		return next_result_idx;
	}

};




/**
 * Sort of external constructor. We want to initialize members of the base class
 * according to the constants in the derived classes.
 *
 * List nodes: Constant, Value have own constructors.
 */
template <class T>
ScalarNode * ScalarExpression::create(ScalarNode *a) {
	T * node_ptr = new T();
	node_ptr->op_code_ = T::op_code;
	node_ptr->add_input(a);
	if (T::n_eval_args == 1) {
		node_ptr->result_storage = none;
	} else {
		ASSERT(T::n_eval_args == 2);
		node_ptr->result_storage = temporary;
	}

	nodes.push_back(node_ptr);
	return node_ptr;
}

template <class T>
ScalarNode * ScalarExpression::create(ScalarNode *a, ScalarNode *b) {
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

	nodes.push_back(node_ptr);
	return node_ptr;
}


} // namespace details
} // namespace bparser
#endif /* INCLUDE_SCALAR_EXPR_HH_ */
