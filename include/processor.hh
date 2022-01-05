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
#include <vector>
#include "config.hh"
#include "assert.hh"
#include "arena_alloc.hh"
#include "expression_dag.hh"
#include "scalar_node.hh"
#include "vectorclass.h"

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


//const uint simd_size = MAX_VECTOR_SIZE / 64;
// const uint simd_size = 4;	
							//myTODO: nesmí být konstanta, ale je třeba nastavit..> v create detekci
							//pote predat do areny jako parametr

//  typedef double double4 __attribute__((__vector_size__(32)));

template <typename VecType>
struct Vec {
	VecType *values;
	uint *subset;

	void set(VecType * v, uint * s) {
		values = v;
		subset = s;
	}

	inline VecType * value(uint i) {
//		std::cout << "self: " << this << std::endl;
//		std::cout << "v: " << values << "s: " << subset << std::endl;
//		std::cout << "i: " << i << "j: " << j << std::endl;
//		std::cout << " si: " << subset[i] << std::endl;
//		std::cout << " v: " << values[subset[i]][j] << "\n";
		return &(values[subset[i]]);
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


//static uint testi = 0;

template<uint NParams, class T, typename VecType>
struct EvalImpl;
//{
//	static inline void eval(Operation op, Workspace &w) {};
//};

template <class T, typename VecType>
struct EvalImpl<1, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w) {
		Vec<VecType> v0 = w.vector[op.arg[0]];
		for(uint i=0; i<w.subset_size; ++i) {
			VecType * v0i = v0.value(i);
			T::eval(*v0i);
		}
	}
};


template <class T, typename VecType>
struct EvalImpl<2, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w) {
		Vec<VecType> v0 = w.vector[op.arg[0]];
		Vec<VecType> v1 = w.vector[op.arg[1]];

		//std::cout << testi++ << " * " << "\n";

		for(uint i=0; i<w.subset_size; ++i) {

			//std::cout << "subset: " << i << std::endl;

			VecType * v0i = v0.value(i);
			VecType * v1i = v1.value(i);
			T::eval(*v0i, *v1i);
		}
	}
};


template <class T, typename VecType>
struct EvalImpl<3, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w) {
		Vec<VecType> v0 = w.vector[op.arg[0]];
		Vec<VecType> v1 = w.vector[op.arg[1]];
		Vec<VecType> v2 = w.vector[op.arg[2]];
//		std::cout << "iv0:" << uint(op.arg[0])
//				<< "iv1:" << uint(op.arg[1])
//				<< "iv2:" << uint(op.arg[2]) << std::endl;

		//std::cout << testi++ << " * " << "\n";
		
		for(uint i=0; i<w.subset_size; ++i) {

			//std::cout << "subset: " << i << std::endl;

			VecType *v0i = v0.value(i);
			VecType *v1i = v1.value(i);
			VecType *v2i = v2.value(i);
			T::eval(*v0i, *v1i, *v2i);
		}
	}
};


template <class T, typename VecType>
struct EvalImpl<4, T, VecType> {
	inline static void eval(Operation op,  Workspace<VecType> &w) {
		Vec<VecType> v0 = w.vector[op.arg[0]];
		Vec<VecType> v1 = w.vector[op.arg[1]];
		Vec<VecType> v2 = w.vector[op.arg[2]];
		Vec<VecType> v3 = w.vector[op.arg[3]];
//		std::cout << "iv0:" << uint(op.arg[0])
//				<< "iv1:" << uint(op.arg[1])
//				<< "iv2:" << uint(op.arg[2]) << std::endl;
		for(uint i=0; i<w.subset_size; ++i) {
			VecType *v0i = v0.value(i);
			VecType *v1i = v1.value(i);
			VecType *v2i = v2.value(i);
			VecType *v3i = v3.value(i);
			T::eval(*v0i, *v1i, *v2i, *v3i);
		}
	}
};



#define CODE(OP_NAME) \
	case (OP_NAME::op_code): operation_eval<OP_NAME>(*op); break

// Note: Internal operations are at most binary, N-ary operations are decomposed into simpler.

struct ProcessorSetup {
	uint vector_size;
	uint n_operations;
	uint n_vectors;
	uint n_constants;
};


//typedef std::conditional<INSTRSET >= 7, Vec4d, Vec2d>::type tmp;
//typedef std::conditional<INSTRSET >= 9, Vec8d, tmp>::type VecType;

template<typename VecType>
struct MyVec {
	typedef VecType Vec;
	uint vector_size;
};

struct ProcessorBase {
	virtual void run() = 0;
	virtual void set_subset(std::vector<uint> const &subset) = 0;

	virtual ~ProcessorBase() {

	}
};

typedef MyVec<double> MyDouble;
typedef MyVec<Vec2d> MyVec2d;
typedef MyVec<Vec4d> MyVec4d;
typedef MyVec<Vec8d> MyVec8d;


/**
 * Store and execute generated "bytecode".
 */

template <typename MV>
struct Processor : ProcessorBase {
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

	typedef typename MV::Vec MVec;
	static const uint simd_size = sizeof(MVec) / sizeof(double);


	static Processor *create(std::vector<ScalarNode *> results, uint vector_size) {
		ExpressionDAG se(results);
			
		return create_processor_(se, vector_size);
	}

	
	static Processor *create_processor_(ExpressionDAG &se, uint vector_size) 
	{
		uint simd_bytes = sizeof(double) * simd_size;
		ExpressionDAG::NodeVec & sorted_nodes = se.sort_nodes();
		//std::cout << "n_nodes: " << sorted_nodes.size() << " n_vec: " << se.n_vectors() << "\n";
		uint memory_est =
				align_size(simd_bytes, sizeof(Processor)) +
				align_size(simd_bytes, sizeof(uint) * vector_size) +
				align_size(simd_bytes, se.temp_end * sizeof(Vec<MVec>)) +
				sizeof(double) * vector_size * (se.temp_end - se.values_end) +
				align_size(simd_bytes, sizeof(MVec) * se.constants_end ) +
				align_size(simd_bytes, sizeof(Operation) * (sorted_nodes.size() + 64) )


				;
		ArenaAlloc arena(simd_bytes, memory_est);

		//uint vec_size = (vector_size / simd_size);
		uint vec_size = simd_size * sizeof(double) * sizeof(double);	//lepsi hodit do promenny?
		return arena.create<Processor>(arena, se, vec_size);
	}

	/**
	 * Do not create processor directly, use the static 'create' method
	 *
	 * vec_size : number of simd blocks (double4).
	 */

	Processor(ArenaAlloc arena, ExpressionDAG &se, uint vec_size)
	: arena_(arena)
	{
		workspace_.vector_size = vec_size;
		workspace_.subset_size = 0;
		workspace_.const_subset = arena_.create_array<uint>(vec_size);
		for(uint i=0; i<vec_size;++i) workspace_.const_subset[i] = 0;
		workspace_.vec_subset = (uint *) arena_.allocate(sizeof(uint) * vec_size);
		//std::cout << "&vec_subset: " << &(workspace_.vec_subset) << "\n";
		//std::cout << "aloc vec_subset: " << workspace_.vec_subset << " size: " << vec_size << "\n";

		//std::cout << "In porcessor.hh: vec_size: " << vec_size << ", simd_size: " << simd_size << std::endl;
		std::cout << "\nIn porcessor.hh: \nvec_size: " << vec_size << "\nsimd_size: " << simd_size << "\nsOfDouble: " << sizeof(double) << "\nsOfMVec: " << sizeof(MVec) << "\nse.temp_end: " << se.temp_end << "\nse.values_end: " << se.values_end << "\nse.constants_end: " << se.constants_end << std::endl;


		workspace_.vector = (Vec<MVec> *) arena_.allocate(sizeof(Vec<MVec>) * se.temp_end);
		MVec * temp_base = (MVec *) arena_.allocate(
				sizeof(double) * vec_size * simd_size * (se.temp_end - se.values_end));
		MVec * const_base = (MVec *) arena_.allocate(
				sizeof(MVec) * se.constants_end);
		for(uint i=0; i< se.constants_end; ++i)
			vec_set(i, const_base + i, workspace_.const_subset);

		uint i_tmp = 0;
		for(uint i=se.values_end; i< se.temp_end; ++i, ++i_tmp)
			vec_set(i, temp_base + i_tmp*vec_size, workspace_.vec_subset);

		// value vectors ... setup when processing the nodes, every value node processed exactly once
		// we need the values pointer from these nodes.

		// TODO: seems that n_temporaries is only for single component of the vector operation
		// however that should be enough as we have separate storage for results so
		// actually we need only 2 temporary vectors
		// need a mean to visualize 'se' graph.
	    auto sorted_nodes = se.sort_nodes();
		uint n_operations = sorted_nodes.size();
		program_ = (Operation *) arena_.allocate(sizeof(Operation) * n_operations);

		/**
		 * TODO separate setup of workspace - no dependence on the order
		 * from composition of operations - top sort but no dep. on result_idx_
		 */
		Operation *op = program_;
		for(auto it=sorted_nodes.rbegin(); it != sorted_nodes.rend(); ++it) {
			//se._print_node(*it);
			ScalarNode * node = *it;
			switch (node->result_storage) {
			case constant: {
				double c_val = *node->get_value();
				MVec * c_ptr = workspace_.vector[node->result_idx_].values;
				c_ptr[0] = c_val;
				break;}
				/*
				for(uint j=0; j<simd_size; ++j)
					c_ptr[0][j] = c_val;
				break;}
				*/
			case value:
				vec_set(node->result_idx_, (MVec *)node->get_value(), workspace_.vec_subset);
				break;
			case temporary:
				*op = make_operation(node);
				++op;
				break;
			case none:
				BP_ASSERT(false);
				//*op = make_operation(node);
				//++op;
				break;
			case expr_result:
				vec_set(node->result_idx_, (MVec *)node->get_value(), workspace_.vec_subset);

				*op = make_operation(node);
				++op;

				//ASSERT(node->n_inputs_ == 1);
				//ScalarNode * prev_node = node->inputs_[0];
				//ASSERT(prev_node->result_storage == temporary);
				//workspace_.vector[prev_node->result_idx_].set((double4 *)node->get_value(), workspace_.vec_subset);
//				std::cout << " ir: " << node->result_idx_ << " a0: "
//						<< workspace_.vector[node->result_idx_].values
//						<< "\n";
				break;
			}
			BP_ASSERT(op < program_ + n_operations);
		}
		op->code = ScalarNode::terminate_op_code;


	}

	void vec_set(uint ivec, MVec * v, uint * s) {
		// std::cout << "Set vec: " << ivec << " ptr: " << &(workspace_.vector[ivec]) << " v: " << v << " s: " << s <<std::endl;
		workspace_.vector[ivec].set(v, s);
	}

	~Processor() {
		arena_.destroy();
	}

	Operation make_operation(ScalarNode * node) {
		Operation op = {(unsigned char)0xff, {0,0,0}}  ;
		op.code = node->op_code_;
		uint i_arg = 0;
		//if (node->result_storage == temporary)
		op.arg[i_arg++] = node->result_idx_;
		for(uint j=0; j<node->n_inputs_; ++j)
			op.arg[i_arg++] = node->inputs_[j]->result_idx_;
		return op;
	}


	template<class T>
	inline void operation_eval(Operation op) {
		EvalImpl<T::n_eval_args, T, MVec>::eval(op, workspace_);
	}

	void run() {
		for(Operation * op = program_;;++op) {
//			std::cout << "op: " << (int)(op->code)
//					<< " ia0: " << (int)(op->arg[0])
//					<< " a0: " << workspace_.vector[op->arg[0]].values
//					<< " ia1: " << (int)(op->arg[1])
//					<< " a1: " << workspace_.vector[op->arg[1]].values
//					<< " ia2: " << (int)(op->arg[2])
//					<< " a2: " << workspace_.vector[op->arg[2]].values << "\n";

			switch (op->code) {
			CODE(_minus_);
			CODE(_add_);
			CODE(_sub_);
			CODE(_mul_);
			CODE(_div_);
			/*
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
			*/
			//CODE(_pow_);
			/*
			CODE(_max_);
			CODE(_min_);
			CODE(_copy_);
			CODE(_ifelse_);
			*/
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
		//std::cout << "vec_subset: " << workspace_.vec_subset << "\n";
		for(uint i=0; i<workspace_.subset_size; ++i) {
			//std::cout << "vec_i: " << workspace_.vec_subset + i << " " << i << "\n";
			workspace_.vec_subset[i] = subset[i];
		}
		// std::cout << "subset: " << workspace_.vec_subset << std::endl;
	}
	

	ArenaAlloc arena_;
	Workspace<MVec> workspace_;
	Operation * program_;
};


} // bparser namespace



#endif /* INCLUDE_PROCESSOR_HH_ */
