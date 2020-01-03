/*
 * processor.hh
 *
 *  Created on: Dec 29, 2019
 *      Author: jb
 */

#ifndef INCLUDE_PROCESSOR_HH_
#define INCLUDE_PROCESSOR_HH_

#include <vector>
#include "config.hh"
#include "vec.hh"
#include "expr.hh"


namespace bparser {


	/**
	 * Memory aligned representation of single operation.
	 */
	struct Operation {
		// Storage pointers
		double * arg_0_ptr;
		double * arg_1_ptr;
		double * result_ptr;
		// Op code. See expr::_xyz_::op_code;
		unsigned char code;
		// Positions of SIMD blocks relative to xyz_ptr.
		unsigned char result_subset;
		unsigned char arg_0_subset;
		unsigned char arg_1_subset;

	};


#define CODE_0_ARY_OP(NAME) \
	case (expr.NAME.op_code) expr.NAME.eval(Vec(op.result_ptr, op.result_subset)); break
#define CODE_1_ARY_OP(NAME) \
	case (expr.NAME.op_code) expr.NAME.eval(Vec(op.result_ptr, op.result_subset), Vec(op.arg_0_ptr, op.arg_0_subset)); break
#define CODE_2_ARY_OP(NAME) \
	case (expr.NAME.op_code) expr.NAME.eval(Vec(op.result_ptr, op.result_subset), Vec(op.arg_0_ptr, op.arg_0_subset), Vec(op.arg_1_ptr, op.arg_1_subset)); break

// Note: Internal operations are at most binary, N-ary operations are decomposed into simpler.


	/**
	 * Store and execute generated "bytecode".
	 */
	struct Processor {
		Processor() {
		}

		void clear(Vec result) {
			program_.clear();
			result = result_;

		}

		void add_op(unsigned char )

		void run() {
			for(const Operation & op : program_) {
				switch (op.code) {
				CODE_0_ARY_OP(_const_);
				CODE_0_ARY_OP(_value_);

				CODE_1_ARY_OP(_abs_);
//				CODE_1_ARY_OP(_exp_);
//				CODE_1_ARY_OP(_pow2_);
//				CODE_1_ARY_OP(_pow10_);
//				CODE_1_ARY_OP(_log_);
//				CODE_1_ARY_OP(_log2_);
//				CODE_1_ARY_OP(_log10_);
//				CODE_1_ARY_OP(_sin_);
//				CODE_1_ARY_OP(_cos_);
//				CODE_1_ARY_OP(_tan_);
//				CODE_1_ARY_OP(_asin_);
//				CODE_1_ARY_OP(_acos_);
//				CODE_1_ARY_OP(_atan_);

				CODE_2_ARY_OP(_add_);
				CODE_2_ARY_OP(_sub_);

				}
			}
		}

		static inline Vec Vec(const double *base_ptr, uint i_subset) {
			return Vec(base_ptr, workspace_.subset(i_subset));
		}

		Workspace workspace_;
		std::vector<Operation> program_;
		Vec result_;
	};

	/**
	 * Create the processor itself and its temporaries in the single chunk of memory.
	 * TODO: Implement and use linear allocator. Can simplify implementation.
	 *
	 * Constraints:
	 * - processor can be applied to different input vectors
	 * - input vectors are at fixed places with fixed maximal size
	 * - the actually used subset of these vectors can change
	 * - the actual size of temporaries change as well
	 * - thus we must reserve the maximal size for the temporaries
	 */
	Processor * create_processor(vector_size) {

	}

} // bparser namespace



#endif /* INCLUDE_PROCESSOR_HH_ */
