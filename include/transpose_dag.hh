/*
 * transpose_dag.hh
 *
 *  Created on: May 5th, 2026
 *      Author: LV
 */

#ifndef INCLUDE_TRANSPOSE_DAG_HH_
#define INCLUDE_TRANSPOSE_DAG_HH_

#include "expression_dag.hh"

namespace bparser {
namespace details {

	struct TransposeNode {
		typedef std::shared_ptr<TransposeNode> TransposeNodePtr;
		typedef std::weak_ptr<TransposeNode> TransposeNodeWPtr;
		ScalarNodePtr node; //node->inputs_
		
		uint n_inputs;
		std::array<TransposeNodePtr, 3> inputs{};
		std::vector<TransposeNodeWPtr> outputs{};
		//weak_ptr should only be null with improper graph handling, since there is a shared_ptr pointing the other way

		TransposeNode(ScalarNodePtr node_ptr) 
		:
			node(node_ptr), 
			n_inputs(node_ptr->n_inputs_)
		{
			/*for (size_t i = 0; i < n_inputs; i++)
			{
				inputs[i] = node->inputs_[i];
			}*/
		}

		TransposeNode(ScalarNodePtr node_ptr, TransposeNodePtr input)
			: TransposeNode(node_ptr)
		{
			inputs[0] = input;
		}

		TransposeNode(ScalarNodePtr node_ptr, TransposeNodePtr input0, TransposeNodePtr input1)
			: TransposeNode(node_ptr)
		{
			inputs[0] = input0;
			inputs[1] = input1;
		}

		TransposeNode(ScalarNodePtr node_ptr, TransposeNodePtr input0, TransposeNodePtr input1, TransposeNodePtr input2)
			: TransposeNode(node_ptr)
		{
			inputs[0] = input0;
			inputs[1] = input1;
			inputs[2] = input2;
		}

		uint n_outputs() const {
			return outputs.size();
		}

		const ScalarNodePtr* inputs_() const{
			return node->inputs_;
		}

		uint n_inputs_() const {
			return node->n_inputs_;
		}

		template<typename T>
		static TransposeNodePtr create() {
			ScalarNodePtr node = ScalarNode::create<T>();
			return std::make_shared<TransposeNode>(node);
		}

		template<typename T>
		static TransposeNodePtr create(TransposeNodePtr input0) {
			ScalarNodePtr node = ScalarNode::create<T>(
				input0->node
				);
			return std::make_shared<TransposeNode>(node,input0);
		}

		template<typename T>
		static TransposeNodePtr create(TransposeNodePtr input0, TransposeNodePtr input1) {
			ScalarNodePtr node = ScalarNode::create<T>(
				input0->node,
				input1->node
			);
			return std::make_shared<TransposeNode>(node,input0, input1);
		}

		template<typename T>
		static TransposeNodePtr create(TransposeNodePtr input0, TransposeNodePtr input1, TransposeNodePtr input2) {
			ScalarNodePtr node = ScalarNode::create<T>(
				input0->node,
				input1->node,
				input2->node
			);
			return std::make_shared<TransposeNode>(node,input0, input1, input2);
		}


	}; //TransposeNode






	class TransposeDAG {
		using TransposeNodePtr = TransposeNode::TransposeNodePtr;
		typedef std::vector<TransposeNodePtr > NodeVec;

		NodeVec nodes;
		std::vector<ScalarNodePtr> results{};

	public:
		TransposeDAG(ExpressionDAG& dag) 
			: TransposeDAG(dag.sort_nodes())
		{
			;
		}

		TransposeDAG(const ExpressionDAG::NodeVec& nodes)
		{
			this->nodes.reserve(nodes.size());
			std::unordered_map<ScalarNodePtr, TransposeNodePtr> ptr_map{};

			for (const ScalarNodePtr& node : nodes) {
				auto tnode = std::make_shared<TransposeNode>(node);
				this->nodes.push_back(tnode);
				ptr_map.emplace(node, tnode);

				if (node->result_storage == ResultStorage::expr_result) {
					results.push_back(node);
				}
			}

			for (TransposeNodePtr& tnode : this->nodes) {
				const ScalarNodePtr& node = tnode->node;

				for (size_t i = 0; i < node->n_inputs_; i++) {
					const ScalarNodePtr& input = node->inputs_[i];

					TransposeNodePtr& tinput = ptr_map[input];
					tnode->inputs[i] = tinput;
					tinput->outputs.push_back(tnode);
				}
			}
		}

		const std::vector<TransposeNodePtr> get_nodes() const {
			return nodes;
		}

		//Result nodes
		const std::vector<ScalarNodePtr> get_results() const {
			return results;
		}

	}; //TransposeDAG

} //details
} //bparser

#endif //INCLUDE_TRANSPOSE_DAG_HH_