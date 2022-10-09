// #include "processor_4.hh"

#include "processor.hh"
#include "arena_alloc.hh"
#include "expression_dag.hh"

// #include "createProcessor.hh"



namespace bparser{

struct Processor_4 : public ProcessorBase {
	static const uint simd_size = sizeof(Vec4d) / sizeof(double);
	/**
	 * Do not create processor directly, use the static 'create' method
	 *
	 * vec_size : number of simd blocks (double4).
	 */
	Processor_4(ArenaAllocPtr arena, ExpressionDAG &se, uint vec_size)
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
		// std::cout << "vec_size: " << vec_size << "\nsimd_size: " << simd_size << "\nsOfVCLVec: " << sizeof(Vec4d) << "\nsOfDouble: " << sizeof(double) << std::endl;
		// std::cout << "se.temp_end: " << se.temp_end << "\nse.values_end: " << se.values_end << "\nse.constants_end: " << se.constants_end << std::endl;


		workspace_.vector = (Vec<Vec4d> *) arena_->allocate(sizeof(Vec<Vec4d>) * se.temp_end);
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
				// c_ptr[0] = c_val;
				// break;}
				
				for(uint j=0; j<simd_size; ++j) {
					// std::cout << "c_val = " << c_val << std::endl;
					// std::cout << "c_ptr[j] = " << c_ptr[j] << std::endl;
					// std::cout << "&c_ptr[j] = " << &c_ptr[j] << std::endl;
					
					c_ptr[j] = c_val;
					// std::cout << "c_ptr[j]after = " << c_ptr[j] << std::endl;
					// std::cout << std::endl;

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
				//*op = make_operation(node);
				//++op;
				break;
			case expr_result:
				for (uint i=0; i < simd_size; i++)
					vec_set(node->result_idx_, (double *)node->get_value(), workspace_.vec_subset);

				*op = make_operation(node);
				++op;
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

	~Processor_4() {
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
		return op;
	}


	template<class T>
	inline void operation_eval(Operation op) {
		EvalImpl<T::n_eval_args, T, Vec4d>::eval(op, workspace_);
	}

	void run() {
		this->copy_inputs();
		for(Operation * op = program_;;++op) {
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
			workspace_.vec_subset[i] = subset[i] * simd_size;									//set po simd_size 0, 4, 8,....
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
	Workspace<Vec4d> workspace_;
	Operation * program_;
	std::vector< std::shared_ptr<ValueCopyNode> > val_copy_nodes_;
};




inline ProcessorBase * create_processor_AVX2(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena) {
    uint simd_bytes = sizeof(double) * simd_size;
    ExpressionDAG::NodeVec & sorted_nodes = se.sort_nodes();
    uint simd_bytes1 = sizeof(Vec4d);
    std::cout << simd_bytes1 << "!=" << simd_bytes << "\n";
    BP_ASSERT(simd_bytes1 == simd_bytes);
    uint vec_size = (vector_size / simd_size);
    uint est = 
            align_size(simd_bytes, sizeof(Processor_4)) +
            align_size(simd_bytes, sizeof(uint) * vector_size) +
            align_size(simd_bytes, se.temp_end * sizeof(Vec<Vec4d>)) +    // temporaries
            align_size(simd_bytes, sizeof(Vec4d) * vec_size * (se.temp_end - se.values_copy_end)) +  // vec_copy, same as temporaries
            align_size(simd_bytes, sizeof(Vec4d) * vec_size * (se.values_copy_end - se.values_end)) + // vector values (probably not neccessary to allocate)
            align_size(simd_bytes, sizeof(Vec4d) * se.constants_end ) +
            align_size(simd_bytes, sizeof(Operation) * (sorted_nodes.size() + 64) );

	// est *= 2;
	// std::cout << "Estimated memory in processor: " << est << std::endl;

    if (arena == nullptr)
        arena = std::make_shared<ArenaAlloc>(simd_bytes, est);
    else
        BP_ASSERT(arena->size_ >= est);
    return arena->create<Processor_4>(arena, se, vec_size);
}

}