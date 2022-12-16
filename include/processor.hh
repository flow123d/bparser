/*
 * processor.hh
 *
 *  Created on: Dec 29, 2019
 *      Author: jb
 */

#ifndef INCLUDE_PROCESSOR_HH_
#define INCLUDE_PROCESSOR_HH_

#include <stdint.h>
#include <malloc.h>
#include <string.h>
#include <vector>
#include "config.hh"
#include "assert.hh"
#include "arena_alloc.hh"
#include "expression_dag.hh"
#include "scalar_node.hh"

namespace bparser {
using namespace details;

/**
 *
 */
//class Workspace {
//public:
//	// Size of the single vector operation. E.g. 4 doubles for AVX2.
//	static const uint simd_block_size = 4;
//
//	Workspace(uint vec_size, uint n_vectors, uint n_constants)
//	: Vec_size_(0)
//	{
//		workspace_size_ = 64;
//		workspace_ = new double[workspace_size_];
//		clear();
//	}
//
//	~Workspace() {
//		delete [] workspace_;
//	}
//
//
//	static void set_subset(std::initializer_list<uint> subset, uint Vec_size) {
//		Workspace::instance().set_subset_(std::vector<uint>(subset), Vec_size);
//	}
//
//	// Set new subset structure and size of the full Vec.
//	static void set_subset(const std::vector<uint> &subset, uint Vec_size) {
//		Workspace::instance().set_subset_(subset, Vec_size);
//	}
//
//	static void set_workspace(uint n_doubles) {
//		Workspace::instance().set_workspace_(n_doubles);
//	}
//
//	// Release all temporary slots.
//	static void clear() {
//		Workspace::instance().next_slot_ = Workspace::instance().workspace_;
//	}
//
//	static uint size() {
//		return Workspace::instance().subset_.size();
//	}
//
//	static double * get_slot()  {
//		return Workspace::instance().get_slot_();
//	}
//
//
//protected:
//	static inline Workspace &instance() {
//		static Workspace w;
//		return w;
//	}
//
//	inline void set_subset_(const std::vector<uint> &subset, uint Vec_size) {
//		clear();
//		Vec_size_ = Vec_size;
//		subset_ = subset; // TODO: avoid copy
//		const_subset_.reserve(size());
//		flat_subset_.reserve(size());
//		for(uint i=0; i<size(); ++i) {
//			const_subset_[i] = 0;
//			flat_subset_[i] = i;
//		}
//	}
//
//	void set_workspace_(uint n_doubles) {
//		delete [] workspace_;
//		workspace_size_ = n_doubles;
//		workspace_ = new double[n_doubles];
//	}
//
//	double * get_slot_()  {
//		ASSERT(next_slot_ < workspace_ + workspace_size_);
//		double *ptr = next_slot_;
//		next_slot_ += size();
//		return ptr;
//	}
//
//
//	static uint * subset() {
//		return &(Workspace::instance().subset_[0]);
//	}
//
//	static uint * flat_subset() {
//		return &(Workspace::instance().flat_subset_[0]);
//	}
//
//	static uint * const_subset() {
//		return &(Workspace::instance().const_subset_[0]);
//	}
//
//
//	std::vector<uint > const_subset_;
//	std::vector<uint > flat_subset_;
//	std::vector<uint > subset_;
//	uint Vec_size_;
//
//	double *workspace_;
//	uint workspace_size_;
//	double *next_slot_;
//};


// const uint simd_size = 4;

// typedef double double4 __attribute__((__vector_size__(32)));

static uint get_simd_size()
{
	if (__builtin_cpu_supports("avx512f"))
	{
		return 8;
	}
	if (__builtin_cpu_supports("avx2"))
	{
		return 4;
	}
	if (__builtin_cpu_supports("sse4.1"))
	{
		return 2;
	}
	else
	{
		return 1;
	}
}

template<typename T>
static T get_true_value();

template<typename T>
static T get_false_value();

template<typename T>
T get_true_value()
{
	T x = 0;
	return as_double(x == x);
}
template<>
double get_true_value<double>()
{
	return double_true();
}

template<typename T>
T get_false_value() {
	T x = 0;
	return as_double(x != x);
}
template<>
double get_false_value<double>()
{
	return double_false();
}

typedef std::shared_ptr<ArenaAlloc> ArenaAllocPtr;

template <typename VecType>
struct Vec {
	double *values;
	uint *subset;

	typedef VecType MyVCLVec;

	void set(double * v, uint * s) {
		values = v;
		subset = s;
	}

	inline double * value(uint i) {
		// std::cout << "self: " << this << std::endl;
		// std::cout << "v: " << values << " s: " << subset << std::endl;
		// std::cout << "i: " << i << std::endl;
		// std::cout << " si: " << subset[i] << std::endl;

		return &(values[subset[i]]);
	}


	VecType true_value() {
		return get_true_value<VecType>();
	}

	VecType false_value() {
		return get_false_value<VecType>();
	}

};


/**
 * Processor's storage.
 */
template <typename VecType>
struct Workspace {
	uint vector_size;

	// Array of vectors. Temporaries, input vectors and result vectors.
	Vec<VecType> *vector;

	uint subset_size;
	uint *const_subset;
	uint *vec_subset;

};


/**
 * Memory aligned representation of single operation.
 */
struct Operation {
	// Op code. See scalar_expr.hh: XYZNode::op_code;
	unsigned char code;
	// index of arguments in the Processors's workspace
	unsigned char arg[4];
};


// static uint testi = 0;

template<uint NParams, class T, typename VecType>
struct EvalImpl;
//{
//	static inline void eval(Operation op, Workspace &w) {};
//};


// EvalmImpl with 1 operand
template <class T, typename VecType>
struct EvalImpl<1, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<1, T, double> {
	inline static void eval(Operation op,  Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<1, T, VecType>::eval(Operation op,  Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		VecType v0i;

		// load value into vector
		v0i.load(v0id);

		// evaluate result
		T::eval(v0i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<1, T, double>::eval(Operation op,  Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		
		// evaluate result
		T::eval(*v0id);
	}
}


// EvalmImpl with 2 operands
template <class T, typename VecType>
struct EvalImpl<2, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<2, T, double> {
	inline static void eval(Operation op,  Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<2, T, VecType>::eval(Operation op,  Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		VecType v0i;
		VecType v1i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);

		// evaluate result
		T::eval(v0i, v1i);

		// store result into memory at v0id
		v0i.store(v0id); 
	}
}

template <class T>
inline void EvalImpl<2, T, double>::eval(Operation op,  Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id);
	}
}


// EvalmImpl with 3 operands
template <class T, typename VecType>
struct EvalImpl<3, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<3, T, double> {
	inline static void eval(Operation op,  Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<3, T, VecType>::eval(Operation op,  Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];
	Vec<VecType> v2 = w.vector[op.arg[2]];
		// std::cout << "iv0:" << uint(op.arg[0])
		// 		<< "iv1:" << uint(op.arg[1])
		// 		<< "iv2:" << uint(op.arg[2]) << std::endl;

	for(uint i=0; i<w.subset_size; ++i) {
		// std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		VecType v0i;
		VecType v1i;
		VecType v2i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);
		v2i.load(v2id);

		// evaluate result
		T::eval(v0i, v1i, v2i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<3, T, double>::eval(Operation op,  Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];
	Vec<double> v2 = w.vector[op.arg[2]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id, *v2id);
	}
}


// EvalmImpl with 3 operands
template <class T, typename VecType>
struct EvalImpl<4, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<4, T, double> {
	inline static void eval(Operation op,  Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<4, T, VecType>::eval(Operation op,  Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];
	Vec<VecType> v2 = w.vector[op.arg[2]];
	Vec<VecType> v3 = w.vector[op.arg[3]];
		// std::cout << "iv0:" << uint(op.arg[0])
		// 		<< "iv1:" << uint(op.arg[1])
		// 		<< "iv2:" << uint(op.arg[2])
		// 		<< "iv3:" << uint(op.arg[3]) << std::endl;

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		double * v3id = v3.value(i);
		VecType v0i;
		VecType v1i;
		VecType v2i;
		VecType v3i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);
		v2i.load(v2id);
		v3i.load(v3id);

		// evaluate result
		T::eval(v0i, v1i, v2i, v3i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<4, T, double>::eval(Operation op,  Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];
	Vec<double> v2 = w.vector[op.arg[2]];
	Vec<double> v3 = w.vector[op.arg[3]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		double * v3id = v3.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id, *v2id, *v3id);
	}
}



#define CODE(OP_NAME) \
	case (OP_NAME::op_code): operation_eval<OP_NAME>(*op); break

// Note: Internal operations are at most binary, N-ary operations are decomposed into simpler.

struct ProcessorSetup {
	uint vector_size;
	uint n_operations;
	uint n_vectors;
	uint n_constants;
};


struct ProcessorBase {
	virtual void run() = 0;
	virtual void set_subset(std::vector<uint> const &subset) = 0;

	ProcessorBase(ArenaAllocPtr arena)
	: arena_(arena) {
		
	}

	virtual ~ProcessorBase() {
	}

	virtual ArenaAllocPtr get_arena(){
		return arena_;
	}
	
	inline static ProcessorBase *create_processor(ExpressionDAG &se, uint vector_size, uint simd_size = 0, ArenaAllocPtr arena = nullptr);

	ArenaAllocPtr arena_;
};



/**
 * Store and execute generated "bytecode".
 */

template <typename VecType>
struct Processor : public ProcessorBase {
	/**
	 *hh
	 * vector_size: maximum vector size in doubles
	 *
	 * TODO: reimplement full__ns to perform topological sort of nodes
	 * TODO: enclose global expression manipulations into a class ScalarExpression
	 * - topological sort
	 * - make_node<> can be its method
	 * - destruction
	 * - extract allocation info
	 * - dependencies
	 * - assigne result ids
	 * - create processor
	 */
	typedef typename VecType::MyVCLVec VCLVec;
	static const uint simd_size = sizeof(VCLVec) / sizeof(double);



// 	static Processor *create(std::vector<ScalarNodePtr > results, uint vector_size) {
// 		ExpressionDAG se(results);
// 			
// 		return create_processor_(se, vector_size);
// 	}



	/**
	 * Do not create processor directly, use the static 'create' method
	 *
	 * vec_size : number of simd blocks (double4).
	 */
	Processor(ArenaAllocPtr arena, ExpressionDAG &se, uint vec_size)
	: ProcessorBase(arena)
	{
		workspace_.vector_size = vec_size;
		workspace_.subset_size = 0;
		workspace_.const_subset = arena_->create_array<uint>(vec_size);
		for(uint i=0; i<vec_size;++i) workspace_.const_subset[i] = 0;
		workspace_.vec_subset = (uint *) arena_->allocate(sizeof(uint) * vec_size);
		
		// std::cout << "&vec_subset: " << &(workspace_.vec_subset) << "\n";
		// std::cout << "aloc vec_subset: " << workspace_.vec_subset << " size: " << vec_size << "\n";

		// std::cout << std::endl << "In porcessor.hh: " << std::endl;
		// std::cout << "vec_size: " << vec_size << "\nsimd_size: " << simd_size << "\nsOfVCLVec: " << sizeof(VCLVec) << "\nsOfDouble: " << sizeof(double) << std::endl;
		// std::cout << "se.temp_end: " << se.temp_end << "\nse.values_end: " << se.values_end << "\nse.constants_end: " << se.constants_end << std::endl;


		workspace_.vector = (Vec<VCLVec> *) arena_->allocate(sizeof(Vec<VCLVec>) * se.temp_end);
		double * temp_base = (double *) arena_->allocate(
				sizeof(double) * vec_size * simd_size * (se.temp_end - se.values_end));
		double * const_base = (double *) arena_->allocate(
				sizeof(double) * simd_size * se.constants_end);
		for(uint i=0; i< se.constants_end; ++i)
			vec_set(i, const_base + i * simd_size, workspace_.const_subset);

		uint i_tmp = 0;
		for(uint i=se.values_end; i< se.values_copy_end; ++i, ++i_tmp)
			vec_set(i, temp_base + i_tmp*vec_size*simd_size, workspace_.vec_subset);

		for(uint i=se.values_copy_end; i< se.temp_end; ++i, ++i_tmp)
			vec_set(i, temp_base + i_tmp*vec_size*simd_size, workspace_.vec_subset);

		// value vectors ... setup when processing the nodes, every value node processed exactly once
		// we need the values pointer from these nodes.

		// TODO: seems that n_temporaries is only for single component of the vector operation
		// however that should be enough as we have separate storage for results so
		// actually we need only 2 temporary vectors
		// need a mean to visualize 'se' graph.
	    auto sorted_nodes = se.sort_nodes();
		uint n_operations = sorted_nodes.size();
		program_ = (Operation *) arena_->allocate(sizeof(Operation) * n_operations);

		/**
		 * TODO separate setup of workspace - no dependence on the order
		 * from composition of operations - top sort but no dep. on result_idx_
		 */
		Operation *op = program_;

		for(auto it=sorted_nodes.rbegin(); it != sorted_nodes.rend(); ++it) {
			// se._print_node(*it);
			// std::cout << "op points at:" << op << std::endl;
			ScalarNodePtr  node = *it;
			switch (node->result_storage) {
			case constant: {
				double c_val = *node->get_value();
				double * c_ptr = workspace_.vector[node->result_idx_].values;
				// std::cout << "node->result_idx_ = " << node->result_idx_ << std::endl;
				// std::cout << "c_ptr = " << c_ptr << std::endl;
				
				for(uint j=0; j<simd_size; ++j) {
					c_ptr[j] = c_val;
				}
				break;}

			case constant_bool:
			{
				double c_val = *node->get_value();
				double * c_ptr = workspace_.vector[node->result_idx_].values;
				Vec<double> v;

				if (c_val == 0.0) {
					for(uint j=0; j<simd_size; ++j)
						c_ptr[j] = v.false_value();
				}
				else {
					for(uint j=0; j<simd_size; ++j)
						c_ptr[j] = v.true_value();
				}
				break;
			}
			case value:
				for (uint i=0; i < simd_size; i++)
					vec_set(node->result_idx_, (double *)node->get_value(), workspace_.vec_subset);
				break;
			case value_copy:
			{
				auto val_copy_ptr = ( std::dynamic_pointer_cast<ValueCopyNode> (node) );
				if (val_copy_ptr->values_ == nullptr) {
					val_copy_ptr->values_ = arena_->create_array<double>(vec_size);
					for (uint i=0; i < simd_size; i++)
						vec_set(node->result_idx_, (double *)node->get_value(), workspace_.vec_subset);
				}
				val_copy_nodes_.push_back(val_copy_ptr);
				break;
			}
			case temporary:
				*op = make_operation(node);
				++op;
				break;
			case none:
				BP_ASSERT(false);
				// *op = make_operation(node);
				// ++op;
				break;
			case expr_result:
				for (uint i=0; i < simd_size; i++)
					vec_set(node->result_idx_, (double *)node->get_value(), workspace_.vec_subset);

				*op = make_operation(node);
				++op;

				// ASSERT(node->n_inputs_ == 1);
				// ScalarNodePtr  prev_node = node->inputs_[0];
				// ASSERT(prev_node->result_storage == temporary);
				// workspace_.vector[prev_node->result_idx_].set((double4 *)node->get_value(), workspace_.vec_subset);
				// std::cout << " ir: " << node->result_idx_ << " a0: "
				// 		<< workspace_.vector[node->result_idx_].values
				// 		<< "\n";
				break;
			}
			BP_ASSERT(op < program_ + n_operations);
		}
		op->code = ScalarNode::terminate_op_code;


	}

	void vec_set(uint ivec, double * v, uint * s) {
		// std::cout << "Set vec: " << ivec << " ptr: " << &(workspace_.vector[ivec]) << " v: " << v  << " &v: " << *v  << " s: " << s << " &s: " << *s <<std::endl;
		workspace_.vector[ivec].set(v, s);
	}

	ArenaAllocPtr get_arena(){
		return arena_;
	}

	~Processor() {
		for (auto node : val_copy_nodes_) {
			node->values_ = nullptr;
		}
		// arena_->destroy();
	}

	Operation make_operation(ScalarNodePtr  node) {
		Operation op = {(unsigned char)0xff, {0,0,0,0}}  ;
		op.code = node->op_code_;
		uint i_arg = 0;
		//if (node->result_storage == temporary)
		op.arg[i_arg++] = node->result_idx_;
		for(uint j=0; j<node->n_inputs_; ++j)
			op.arg[i_arg++] = node->inputs_[j]->result_idx_;

		// std::cout << "Created new op: " << (int)(op.code)
		// 	<< " ia0: " << (int)(op.arg[0])
		// 	<< " a0: " << workspace_.vector[op.arg[0]].values
		// 	<< " ia1: " << (int)(op.arg[1])
		// 	<< " a1: " << workspace_.vector[op.arg[1]].values
		// 	<< " ia2: " << (int)(op.arg[2])
		// 	<< " a2: " << workspace_.vector[op.arg[2]].values
		// 	<< " ia3: " << (int)(op.arg[3])
		// 	<< " a3: " << workspace_.vector[op.arg[3]].values << "\n";
		return op;
	}


	template<class T>
	inline void operation_eval(Operation op) {
		EvalImpl<T::n_eval_args, T, VCLVec>::eval(op, workspace_);
	}

	void run() {
		this->copy_inputs();
		for(Operation * op = program_;;++op) {
			// std::cout << "op points at:" << op << std::endl;

			// std::cout << "op: " << (int)(op->code)
			// 		<< " ia0: " << (int)(op->arg[0])
			// 		<< " a0: " << workspace_.vector[op->arg[0]].values
			// 		<< " ia1: " << (int)(op->arg[1])
			// 		<< " a1: " << workspace_.vector[op->arg[1]].values
			// 		<< " ia2: " << (int)(op->arg[2])
			// 		<< " a2: " << workspace_.vector[op->arg[2]].values
			// 		<< " ia3: " << (int)(op->arg[3])
			// 		<< " a3: " << workspace_.vector[op->arg[3]].values << "\n";

			switch (op->code) {
			CODE(_minus_);
			CODE(_add_);
			CODE(_sub_);
			CODE(_mul_);
			CODE(_div_);
			CODE(_mod_);
			CODE(_eq_);
			CODE(_ne_);
			CODE(_lt_);
			CODE(_le_);
			CODE(_neg_);
			CODE(_or_);
			CODE(_and_);
			CODE(_abs_);
			CODE(_sqrt_);
			CODE(_exp_);
			CODE(_log_);
			CODE(_log10_);
			CODE(_sin_);
			CODE(_sinh_);
			CODE(_asin_);
			CODE(_cos_);
			CODE(_cosh_);
			CODE(_acos_);
			CODE(_tan_);
			CODE(_tanh_);
			CODE(_atan_);
			CODE(_ceil_);
			CODE(_floor_);
			CODE(_isnan_);
			CODE(_isinf_);
			CODE(_sgn_);
			CODE(_atan2_);
			CODE(_pow_);
			CODE(_max_);
			CODE(_min_);
			CODE(_copy_);
			CODE(_ifelse_);
			CODE(_log2_);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
//			CODE(__);
			case (ScalarNode::terminate_op_code): return; // terminal operation
			}
		}
	}

	// Set subset indices of active double4 blocks.
	// TODO: Provide getter for pointer to the workspace subset in order to
	// fill it (some where), can be passed together with fixed size as std::span
	void set_subset(std::vector<uint> const &subset)
	{
		BP_ASSERT( (subset.size() <= workspace_.vector_size) );
		workspace_.subset_size = subset.size();
		// std::cout << "vec_subset: " << workspace_.vec_subset << "\n";
		for(uint i=0; i<workspace_.subset_size; ++i) {
			// std::cout << "subset_i: " << subset[i] << " i=" << i << "\n";
			workspace_.vec_subset[i] = subset[i] * simd_size;
			// std::cout << "subsetvec_i: " << workspace_.vec_subset[i]<< " i=" << i << "\n";
		}
		// std::cout << "subset: " << workspace_.vec_subset << std::endl;
	}

	// Copy data of ValueCopyNode objects to arena_
	void copy_inputs()
	{

		for (auto node : val_copy_nodes_) {
			memcpy(node->values_, node->source_ptr_, workspace_.vector_size * sizeof *node->values_);
		}
	}

	// ArenaAlloc arena_;
	Workspace<VCLVec> workspace_;
	Operation * program_;
	std::vector< std::shared_ptr<ValueCopyNode> > val_copy_nodes_;
};


template <class VCLVec> 
ProcessorBase * create_processor_(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena)
{
    uint simd_bytes = sizeof(double) * simd_size;
    ExpressionDAG::NodeVec & sorted_nodes = se.sort_nodes();
    uint simd_bytes1 = sizeof(VCLVec);
    // std::cout << simd_bytes1 << "!=" << simd_bytes << "\n";
    BP_ASSERT(simd_bytes1 == simd_bytes);
    uint vec_size = (vector_size / simd_size);
    uint est = 
            align_size(simd_bytes, sizeof(Processor<Vec<VCLVec>>)) +
            align_size(simd_bytes, sizeof(uint) * vector_size) +
            align_size(simd_bytes, se.temp_end * sizeof(Vec<VCLVec>)) +
            align_size(simd_bytes, sizeof(VCLVec) * vec_size * (se.temp_end - se.values_copy_end)) +  // vec_copy, same as temporaries
            align_size(simd_bytes, sizeof(VCLVec) * vec_size * (se.values_copy_end - se.values_end)) + // vector values (probably not neccessary to allocate)
            align_size(simd_bytes, sizeof(VCLVec) * se.constants_end ) +
            align_size(simd_bytes, sizeof(Operation) * (sorted_nodes.size() + 64) );

	// std::cout << "Estimated memory in processor: " << est << std::endl;

    if (arena == nullptr)
        arena = std::make_shared<ArenaAlloc>(simd_bytes, est);
    else
        BP_ASSERT(arena->size_ >= est);
    return arena->create<Processor<Vec<VCLVec>>>(arena, se, vec_size);
}


} // bparser namespace


#endif /* INCLUDE_PROCESSOR_HH_ */
