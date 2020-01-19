/*
 * expression_dag.hh
 *
 *  Created on: Jan 19, 2020
 *      Author: jb
 */

#ifndef INCLUDE_EXPRESSION_DAG_HH_
#define INCLUDE_EXPRESSION_DAG_HH_


#include <vector>
#include <cmath>
#include <map>
#include "config.hh"
#include "assert.hh"
#include "scalar_expr.hh"



namespace bparser {
namespace details {




/**
 * Auxiliary class to form evaluation graph (DAG) and make topological sort
 * to get order of operations.
 * TODO: optimize topological sort for the number of temporaries.
 */
class ExpressionDAG {
public:
	typedef std::vector<ScalarNode *> NodeVec;

private:
	/// All nodes in the expressions (reached from the results).
	NodeVec nodes;
	/// Backward topologicaly sorted nodes; results first, inputs last
	NodeVec sorted;
	/// Result nodes, given as input.
	NodeVec results;


	/**
	 * Used in the setup_result_storage to note number of unclosed nodes
	 * dependent on the value. When this drops to zero the temporary may be reused.
	 * TODO: use to reorder nodes in the topological sort to minimize number of temporaries
	 */
	std::vector<uint> storage;

public:
	/**
	 * All input, temporary and result vectors (represented by bparser::details::Vec)
	 * are assigned to the fixed position in the storage table of the Workspace.
	 * First are constants, then external vectors (Value and Result) finally temporaries.
	 */
	/// End index of the constant vectors in storage.
	uint constants_end;
	/// End index of the values and result vectors in the storage.
	uint values_end;
	/// End index of the temporary vectors in storage.
	uint temp_end;


	ExpressionDAG(std::vector<ScalarNode *> res)
	:
		results(res.begin(), res.end()),
		constants_end(0),
		values_end(0),
		temp_end(0)
	{
		sort_nodes();
	}


	~ExpressionDAG() {
		for(ScalarNode *node : nodes) {
			delete node;
		}
	}

	/**
	 * Return nodes in the topological order (result nodes first).
	 * It also assign position of the node results in the storage (result_idx_).
	 */
	NodeVec & sort_nodes() {
        if (sorted.size() > 0)
        	return sorted;

        /**
         * TODO: there is some infinite loop
         */
        _collect_nodes();
		ASSERT(sorted.size() == 0);
		_topological_sort();

		_setup_result_storage();
		temp_end += storage.size();
		return sorted;
	}


	/**
	 * Print ScalarExpression graph in the dot format.
	 */
	void print_in_dot() {
		std::map<ScalarNode *, uint> i_node;
		sort_nodes();
		for(uint i=0; i<sorted.size(); ++i) i_node[sorted[i]] = i;


		std::cout << "\n" << "----- begin cut here -----" << "\n";
		std::cout << "digraph Expr {" << "\n";
		for(uint i=0; i<sorted.size(); ++i) {
			i_node[sorted[i]] = i;
			//std::cout << i << " n: " << sorted[i]->n_inputs_ << "\n";
			for(uint in=0; in<sorted[i]->n_inputs_; ++in ) {
				std::cout << "    ";
				_print_node(sorted[i]);
				std::cout << " -> ";
				_print_node(sorted[i]->inputs_[in]);
				std::cout << "\n";
			}
		}
		std::cout << "}" << "\n";
		std::cout << "\n" << "----- end   cut here -----" << "\n";
		std::cout.flush();
	}

	void _print_node(ScalarNode * node) {
		std::cout << node->op_name_ <<  "_" << node->result_idx_;
	}


private:
	void _print_i_node(uint i) {
		std::cout << sorted[i]->op_name_ << "_" << i << "_"<< sorted[i]->result_idx_;
	}



	/**
	 * Performs BFS to:
	 * - collect all nodes in the expression graph
	 * - count constants and values/results
	 * - assign constant and value result_idx_ to the nodes.
	 */
	void _collect_nodes() {
		// ScalarNode::reslut_idx_ == -1,
		// we set it to -2 to identify passed nodes
		nodes.clear();

		// collect nodes
		for(auto node : results)
			node->result_idx_ = -2;
		nodes.insert(nodes.begin(), results.begin(), results.end());
		for(uint i=0; i < nodes.size(); ++i) {
			ScalarNode * node = nodes[i];
			for(uint in=0; in < node->n_inputs_; ++in)  {
				ScalarNode * other = node->inputs_[in];
				if (other->result_idx_ != -2) {
					ASSERT(other->result_idx_ == -1);
					ASSERT(other->result_storage != expr_result);
					other->result_idx_ = -2;
					nodes.push_back(other);
				}
			}
		}


		// set result_idx_ of constant nodes
		uint i_storage = 0;
		for(ScalarNode * node : nodes)
			if (node->result_storage == constant)
				node->result_idx_ = i_storage++;
		constants_end = i_storage;
		// set result_idx_ of value/result nodes
		for(ScalarNode * node : nodes)
			if (node->result_storage == value || node->result_storage == expr_result)
				node->result_idx_ = i_storage++;
		values_end = i_storage;
		temp_end = i_storage; // still empty
	}


	void _topological_sort() {

		// in-degree of nodes (number of dependent nodes).
		for(ScalarNode *node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNode *node : nodes) {
			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ += 1;
			}
		}

		// Kahn's algorithm for topo. sort
		// Drawing nodes form the stack in different order leads to all possible topological orderings.
		// However stack seems to reuse temporaries more efficiently then e.g. queue.
		// Yet it is not optimal.
		// TODO: probable optimal algorithm viz. TGH semestralky 2020
		NodeVec stack;
		for(auto node: nodes)
			if (node->n_dep_nodes_ == 0)
				stack.push_back(node);

		while (stack.size() > 0) {
			ScalarNode * node = stack.back();
			// std::cout << "node: " << node << " res: " << node->result_storage << " n_dep: " << node->n_dep_nodes_  << "\n";
			stack.pop_back();
			sorted.push_back(node);

			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ -= 1;
				//std::cout << "  node: " << node->inputs_[in] << " n_dep: " << node->inputs_[in]->n_dep_nodes_ << "\n";
				if (node->inputs_[in]->n_dep_nodes_ == 0)
					stack.push_back(node->inputs_[in]);
			}
		}

	}

	/**
	 * Assign result_idx_ to the temporary nodes, reusing
	 * storage positions.
	 */
	void _setup_result_storage() {
		// in-degree of nodes (number of dependent nodes).
		for(ScalarNode *node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNode *node : nodes)
			for(uint in=0; in < node->n_inputs_; ++in) node->inputs_[in]->n_dep_nodes_ += 1;

		// Mimic expression evaluation, reversed topological order.
		for(auto it=sorted.rbegin(); it != sorted.rend(); ++it) {
			ScalarNode * node = *it;
			_allocate_storage(node);
			for(uint in=0; in < node->n_inputs_; ++in) {
				node->inputs_[in]->n_dep_nodes_ -= 1;
				if (node->inputs_[in]->n_dep_nodes_ == 0) {
					_deallocate_storage(node->inputs_[in]);
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
	void _allocate_storage(ScalarNode *node) {
		if (node->result_storage == temporary) {
			for(uint i=0; i<storage.size(); ++i)
				if (storage[i] == 0) {
					storage[i] = 1;
					node->result_idx_ = temp_end + i;
					return;
				}
			node->result_idx_ = temp_end + storage.size();
			storage.push_back(1);
			return;
		}
	}


	void _deallocate_storage(ScalarNode *node) {
		if (node->result_storage == temporary) {
			storage[node->result_idx_ - temp_end] = 0;
		}
	}










};






} // namespace details
} // namespace bparser




#endif /* INCLUDE_EXPRESSION_DAG_HH_ */
