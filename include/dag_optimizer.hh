/*
 * dag_optimizer.hh
 *
 *  Created on: May 5th, 2026
 *      Author: LV
 */

#ifndef INCLUDE_DAG_OPTIMIZER_HH_
#define INCLUDE_DAG_OPTIMIZER_HH_

#include "transpose_dag.hh"

namespace bparser {
namespace details {

	using TransposeNodePtr = TransposeNode::TransposeNodePtr;

	struct DAGOptimization {
		virtual bool can_optimize(TransposeNodePtr) = 0;
		virtual void apply(TransposeNodePtr) = 0;

		virtual ~DAGOptimization() = default;

		template<typename T>
		bool is_op(const TransposeNodePtr& tnode) const {
			return is_op<T>(tnode->node);
		}

		template<typename T>
		bool is_op(const ScalarNodePtr& node) const{
			return node->op_code_ == T::op_code;
		}


		//  current <  - - - a		 current          a
		//			  \			  	               /
		//			    \		 =>		         /
		//  new_node	  \ b	     new_node < - - - b
		//

		void replace_outputs_input(TransposeNodePtr current, ScalarNodePtr new_node) {
			for (TransposeNodePtr output : current->outputs) {
				for (size_t i = 0; i < output->n_inputs_(); i++) {
					if (output->inputs_()[i] == current->node) {
						output->node->inputs_[i] = new_node;
					}
				}
			}
		}
		void replace_outputs_input(TransposeNodePtr current, TransposeNodePtr new_node) {
			replace_outputs_input(current, new_node->node);
		}

		//  a - - - > current         a         current
		//		  /			  	        \\		.
		//		/			  =>		  \\    .  
		//  b /		new_node	      b - - - > new_node
		//

		void replace_inputs_output(TransposeNodePtr current, TransposeNodePtr new_node) {
			for (uint i = 0; i < current->n_inputs; i++) {
				TransposeNodePtr& input = current->inputs[i];
				for (size_t j = 0; j < input->n_outputs(); j++) {
					if (input->outputs[j] == current) {
						input->outputs[j] = new_node;
					}
				}
			}
		}

	};
	typedef std::shared_ptr<DAGOptimization> OptPtr;

	class DAGOptimizer {

		std::vector<OptPtr> optimizations{};
	public:
		DAGOptimizer(const std::vector<OptPtr>& opts) 
			:
			optimizations(opts) 
		{
			;
		}

		//Applies in-place optimizations to the DAG, which invalidates the ExpressionDAG's nodes and sorted fields. Use the returned ExpressionDAG instead!
		ExpressionDAG optimize(ExpressionDAG& dag) {
			TransposeDAG tdag(dag);
			return optimize(tdag);
		}

		//Applies in-place optimizations to the DAG, which invalidates the ExpressionDAG's nodes and sorted fields. Use the returned ExpressionDAG instead!
		ExpressionDAG optimize(TransposeDAG& tdag) {
			std::vector<std::pair<TransposeNodePtr, OptPtr>> possible_opts{};
			for (TransposeNodePtr tnode : tdag.get_nodes()) {
				for (OptPtr& opt : optimizations) {
					if (opt->can_optimize(tnode)) {
						possible_opts.emplace_back(tnode, opt);
					}
				}
			}

			for (auto& [tnode, opt] : possible_opts) {
				if (opt->can_optimize(tnode)) opt->apply(tnode);
			}

			//Recreate the ExpressioDAG, since its sorted vector is no longer valid
			return ExpressionDAG(tdag.get_results());
		}

	};

	//(a * b) + c
	struct MulAddOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return	
				tnode->n_outputs() == 1 &&
				is_op<_mul_>(tnode) &&
				is_op<_add_>(tnode->outputs[0])
				;
		}

		void apply(TransposeNodePtr mul) override {
			TransposeNodePtr add = mul->outputs[0];

			TransposeNodePtr muladd = TransposeNode::create<_muladd_>(
				mul->inputs[0],
				mul->inputs[1],
				add->inputs_()[0] == mul->node ? add->inputs[1] : add->inputs[0]
			);
			muladd->outputs = add->outputs;

			replace_inputs_output(mul, muladd);
			replace_inputs_output(add, muladd);

			replace_outputs_input(add, muladd);
		}
	};

	//(a * b) - c
	struct MulSubOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return
				tnode->n_outputs() == 1 &&
				is_op<_mul_>(tnode) &&
				is_op<_sub_>(tnode->outputs[0]) &&
				tnode->outputs[0]->inputs_()[0] == tnode->node // a*b - c not c - a*b
				;
		}

		void apply(TransposeNodePtr mul) override {
			TransposeNodePtr sub = mul->outputs[0];

			TransposeNodePtr mulsub = TransposeNode::create<_mulsub_>(
				mul->inputs[0],
				mul->inputs[1],
				sub->inputs[1]
			);
			mulsub->outputs = sub->outputs;

			replace_inputs_output(mul, mulsub);
			replace_inputs_output(sub, mulsub);

			replace_outputs_input(sub, mulsub);
		}
	};

} //details
} //bparser

#endif //INCLUDE_DAG_OPTIMIZER_HH_