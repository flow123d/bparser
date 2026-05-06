/*
 * dag_optimizer.hh
 *
 *  Created on: May 3, 2026
 *      Author: LV
 */

#ifndef INCLUDE_DAG_OPTIMIZER_HH_
#define INCLUDE_DAG_OPTIMIZER_HH_

#include "transpose_dag.hh"

namespace bparser {
namespace details {

	using TransposeNodePtr = TransposeNode::TransposeNodePtr;
	using TransposeNodeWPtr = TransposeNode::TransposeNodeWPtr;

	struct DAGOptimization {
		virtual bool can_optimize(TransposeNodePtr) = 0;
		virtual void apply(TransposeNodePtr) = 0;

		virtual ~DAGOptimization() = default;

		template<typename T>
		bool is_op(const TransposeNodePtr& tnode) const {
			return is_op<T>(tnode->node);
		}

		template<typename T>
		bool is_op(const TransposeNodeWPtr& tnode) const {
			TransposeNodePtr ptr = tnode.lock();
			if (!ptr) return false;
			return is_op<T>(ptr->node);
		}

		template<typename T>
		bool is_op(const ScalarNodePtr& node) const{
			return node->op_code_ == T::op_code;
		}

		bool is_result(const TransposeNodePtr& tnode) const {
			return is_result(tnode->node);
		}

		bool is_result(const TransposeNodeWPtr& tnode) const {
			TransposeNodePtr ptr = tnode.lock();
			if (!ptr) return false;
			return is_result(ptr->node);
		}

		bool is_result(const ScalarNodePtr& node) const {
			return node->result_storage == ResultStorage::expr_result;
		}


		//  current <  - - - a		 current          a
		//			  \			  	               /
		//			    \		 =>		         /
		//  new_node	  \ b	     new_node < - - - b
		//

		void replace_outputs_input(TransposeNodePtr current, TransposeNodePtr new_node) {
			for (TransposeNodeWPtr woutput : current->outputs) {
				TransposeNodePtr output = woutput.lock();
				if (!output) continue;
				for (size_t i = 0; i < output->n_inputs_(); i++) {
					if (output->inputs_()[i] == current->node) {
						output->inputs[i] = new_node;
						output->node->inputs_[i] = new_node->node;
					}
				}
			}
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
					if (input->outputs[j].lock()->node == current->node) {
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




	//OPTIMIZATIONS ----------------------------------------------------------------


	//(a * b) + c
	struct MulAddOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return	
				tnode->n_outputs() == 1 &&
				is_op<_mul_>(tnode) &&
				is_op<_add_>(tnode->outputs[0]) &&
				!is_result(tnode->outputs[0])
				;
		}

		void apply(TransposeNodePtr mul) override {
			TransposeNodePtr add = mul->outputs[0].lock();

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
				!is_result(tnode->outputs[0]) &&
				tnode->outputs[0].lock()->inputs_()[0] == tnode->node // a*b - c not c - a*b
				;
		}

		void apply(TransposeNodePtr mul) override {
			TransposeNodePtr sub = mul->outputs[0].lock();

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

	//c - (a * b)
	struct NMulAddOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return
				tnode->n_outputs() == 1 &&
				is_op<_mul_>(tnode) &&
				is_op<_sub_>(tnode->outputs[0]) &&
				!is_result(tnode->outputs[0]) &&
				tnode->outputs[0].lock()->inputs_()[1] == tnode->node //  c - a*b not a*b - c
				;
		}

		void apply(TransposeNodePtr mul) override {
			TransposeNodePtr sub = mul->outputs[0].lock();

			TransposeNodePtr mulsub = TransposeNode::create<_nmuladd_>(
				mul->inputs[0],
				mul->inputs[1],
				sub->inputs[0]
			);	
			mulsub->outputs = sub->outputs;

			replace_inputs_output(mul, mulsub);
			replace_inputs_output(sub, mulsub);

			replace_outputs_input(sub, mulsub);
		}
	};

	//(a + b) * c
	struct AddMulOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return
				tnode->n_outputs() == 1 &&
				is_op<_add_>(tnode) &&
				is_op<_mul_>(tnode->outputs[0]) &&
				!is_result(tnode->outputs[0])
				;
		}

		void apply(TransposeNodePtr add) override {
			TransposeNodePtr mul = add->outputs[0].lock();

			TransposeNodePtr addmul = TransposeNode::create<_addmul_>(
				add->inputs[0],
				add->inputs[1],
				mul->inputs_()[0] == add->node ? mul->inputs[1] : mul->inputs[0]
			);
			addmul->outputs = mul->outputs;

			replace_inputs_output(add, addmul);
			replace_inputs_output(mul, addmul);

			replace_outputs_input(mul, addmul);
		}
	};

	//(a - b) * c
	struct SubMulOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return
				tnode->n_outputs() == 1 &&
				is_op<_sub_>(tnode) &&
				is_op<_mul_>(tnode->outputs[0]) &&
				!is_result(tnode->outputs[0])
				;
		}

		void apply(TransposeNodePtr sub) override {
			TransposeNodePtr mul = sub->outputs[0].lock();

			TransposeNodePtr mulsub = TransposeNode::create<_submul_>(
				sub->inputs[0],
				sub->inputs[1],
				mul->inputs_()[0] == sub->node ? mul->inputs[1] : mul->inputs[0]
			);
			mulsub->outputs = mul->outputs;

			replace_inputs_output(sub, mulsub);
			replace_inputs_output(mul, mulsub);

			replace_outputs_input(mul, mulsub);
		}
	};

	//(a * b) * c
	struct MulMulOpt : public DAGOptimization {
		bool can_optimize(TransposeNodePtr tnode) override {
			return
				tnode->n_outputs() == 1 &&
				is_op<_mul_>(tnode) &&
				is_op<_mul_>(tnode->outputs[0]) &&
				!is_result(tnode->outputs[0])
				;
		}

		void apply(TransposeNodePtr mul0) override {
			TransposeNodePtr mul1 = mul0->outputs[0].lock();

			TransposeNodePtr addmul = TransposeNode::create<_mulmul_>(
				mul0->inputs[0],
				mul0->inputs[1],
				mul1->inputs_()[0] == mul0->node ? mul1->inputs[1] : mul1->inputs[0]
			);
			addmul->outputs = mul1->outputs;

			replace_inputs_output(mul0, addmul);
			replace_inputs_output(mul1, addmul);

			replace_outputs_input(mul1, addmul);
		}
	};

} //details
} //bparser

#endif //INCLUDE_DAG_OPTIMIZER_HH_