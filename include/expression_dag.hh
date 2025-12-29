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
#include "scalar_node.hh"
#include "assert.hh"
#include "array.hh"


namespace bparser {
namespace details {




/**
 * Auxiliary class to form evaluation graph (DAG) and make topological sort
 * to get order of operations.
 * TODO: optimize topological sort for the number of temporaries.
 */
class ExpressionDAG {
public:
	typedef std::vector<ScalarNodePtr > NodeVec;

private:
	/// All nodes in the expressions (reached from the results).
	NodeVec nodes;
	/// Backward topologicaly sorted nodes; results first, inputs last
	NodeVec sorted;
	/// Result nodes, given as input.
	NodeVec results;

	typedef std::pair<std::string, bool> InvDotNameAndScalar;
	typedef std::map<ScalarNodePtr, InvDotNameAndScalar> InvDotMap;
	typedef std::unordered_map<double*, std::string> CXXVarMap;

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
	/// End index of the copied values and result vectors in the storage.
	uint values_copy_end;
	/// End index of the temporary vectors in storage.
	uint temp_end;


	ExpressionDAG(std::vector<ScalarNodePtr > res)
	:
		results(res.begin(), res.end()),
		constants_end(0),
		values_end(0),
		values_copy_end(0),
		temp_end(0)
	{
		sort_nodes();
	}


	~ExpressionDAG() {
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
		BP_ASSERT(sorted.size() == 0);
		_topological_sort();

		_setup_result_storage();
		temp_end += storage.size();
		return sorted;
	}


	/**
	 * Print ScalarExpression graph in the dot format.
	 * Useful for debugging
	 */
	void print_in_dot() {
		std::map<ScalarNodePtr , uint> i_node;
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

	void _print_node(ScalarNodePtr  node) {
		std::cout << "Node: " << node->op_name_ <<  "_" << node->result_idx_ << " " << node->result_storage << std::endl;
	}

	/**
	 * Print ScalarExpression graph in the common dot format.
	 * Useful for understanding the DAG.
	 */
	void print_in_dot2() {
		print_in_dot2(InvDotMap());
	}

	/**
	 * Print ScalarExpression graph in the common dot format.
	 * Useful for understanding the DAG. Using the parser's map of var. Name -> Array find the inverse ScalarNodePtr -> var. Name
	 */
	void print_in_dot2(const std::map<std::string, bparser::Array>& symbols) {
		print_in_dot2(create_inverse_map(symbols));
	}

	/**
	 * Print ScalarExpression graph in the common dot format.
	 * Useful for understanding the DAG. Using the map of ScalarNodePtr -> variableName
	 */
	void print_in_dot2(const InvDotMap& names) {

		sort_nodes();
		
		std::cout << "\n" << "----- begin cut here -----" << "\n";
		std::cout << "digraph Expr {" << "\n";

		std::cout << "/* definitions */" << "\n";

		std::cout << "edge [dir=back]" << "\n";
		for (uint i = 0; i < sorted.size(); ++i) {
			_print_dot_node_definition(sorted[i],names);
		}
		std::cout << "/* end of definitions */" << "\n";

		for (uint i = 0; i < sorted.size(); ++i) {
			for (uint in = 0; in < sorted[i]->n_inputs_; ++in) {
				std::cout << "    ";
				std::cout << _get_dot_node_id(sorted[i]);
				std::cout << "\n -> ";
				std::cout << _get_dot_node_id(sorted[i]->inputs_[in]);
				std::cout << "\n\n";
			}
		}
		std::cout << "}" << "\n";
		std::cout << "-----  end cut here  -----" << "\n";
		std::cout.flush();
	}

	std::string print_in_cxx(const CXXVarMap& map) {
		std::ostringstream result;
		NodeVec result_nodes;

		result << "//AUTOGENERATED This file has been autogenerated by bparser::ExpressionDAG::print_in_cxx" << "\n";
		result << "\n";
		result << "#ifndef NITPICK_INCLUDE_ONLY" << "\n";
		result << "#include \"parser.hh\"" << "\n";
		result << "using namespace bparser;" << "\n";
		result << "using namespace bparser::details;" << "\n";
		result << "int main(){ //This is here only to stop any IDE warnings, do not run this file as is" << "\n";
		result << "#endif //NITPICK_INCLUDE_ONLY" << "\n";
		result << "\n";

		//Print nodes
		for (uint i = sorted.size(); i-- > 0U; ) { // N-1,N-2,... 0
			
			result << _get_cxx_node_definition(sorted[i],map);
			//result << "\n";
			if (sorted[i]->result_storage == expr_result) {
				result_nodes.push_back(sorted[i]);
			}
		}
		result << "\n\n";
		//Print results
		for (uint i = 0U; i < result_nodes.size(); ++i) {
			result << "ScalarNodePtr " << "r" << i << " = " << _get_cxx_result(result_nodes[i],map);
		}
		//Print final dag
		result << "ExpressionDAG se({";
		for (uint i = 0U; i < result_nodes.size(); ++i) {
			result << "r" << i;
			if (i < result_nodes.size() - 1) {
				result << ", ";
			}
		}
		result << "}); //Do not rename this DAG, it is used later in the nitpick_run.cc file";

		result << "\n";
		result << "#ifndef NITPICK_INCLUDE_ONLY" << "\n";
		result << "} //main" << "\n";
		result << "#endif //NITPICK_INCLUDE_ONLY" << "\n";
		//std::cout << result.str();
		return result.str();
		
	}
	
	//Create a map of ScalarNodePtr -> (variable name, is_scalar)
	InvDotMap create_inverse_map(const std::map<std::string, bparser::Array>& symbols) const {
		InvDotMap inv_map;
		if (symbols.empty()) return inv_map;
		for (const auto& s : symbols)
		{
			for (const auto& n : s.second.elements()) {
				inv_map[n] = std::pair<std::string,bool>(s.first, s.second.shape().empty());
			}
		}
		return inv_map;
	}


private:
	//Print the vertice identifier for dot
	std::string _get_dot_node_id(const ScalarNodePtr& node) const {
		return (std::ostringstream() << node->op_name_ << "_" << (uintptr_t)node.get() << "__" << node->result_storage).str();// << std::endl;
	}

	//Print how the vertice should look in dot
	void _print_dot_node_definition(const ScalarNodePtr& node, const InvDotMap& invmap) const {
		std::cout << _get_dot_node_id(node);
		std::cout << ' ';

		if (node->result_storage == ResultStorage::constant) {				// Constant
			std::cout << "[shape=circle,";

			try { //If the constant has a name
				std::string name(invmap.at(node).first);
				std::cout << "label=\"" << name << ": " << *node->values_ << "\",group=\"" << name << '"';
			}
			catch (const std::out_of_range&) { //No name
				std::cout << "label=\"" << "const " << *node->values_ << '"';
			}
			std::cout  << "]" << std::endl;
		}

		else if (node->result_storage == ResultStorage::constant_bool) {	//Constant bool
			std::cout << "[shape=circle,";

			try { //If the constant has a name
				std::string name(invmap.at(node).first);
				std::cout << "label=\"" << name << ": " << *node->values_ << "\",group=\"" << name << '"';
			}
			catch (const std::out_of_range&) { //No name
				std::cout << "label=\"" << "const " << *node->values_ << '"';
			}
			std::cout << "]" << std::endl;
		}

		else if (node->result_storage == ResultStorage::expr_result) {		//Result
			std::cout << "[shape=box,label=\"" << node->op_name_ << " [" << node->result_idx_ << "]" << "\"]" << std::endl;
		}

		else if (node->result_storage == ResultStorage::value) {			// Value

			std::cout << "[shape=circle,";
			try {
				std::string name(invmap.at(node).first);
				bool scalar(invmap.at(node).second);
				if (scalar) {
					std::cout << "label=\"" << name << '"';
				}
				else {
					std::cout << "label=<" << name << "<SUB>i</SUB>" << '>';
				}
				std::cout << ",group=\"" << name << '"';
			}
			catch (const std::out_of_range&) {
				std::cout << "label=<<I>var</I>>";
			}
			
			std::cout << "]" << std::endl;
		}

		else if (node->result_storage == ResultStorage::value_copy) {		//Value copy
			std::cout << "[shape=circle,";
			try {
				std::string name(invmap.at(node).first);
				bool scalar(invmap.at(node).second);
				if (scalar) {
					std::cout << "label=\"" << name << '"';
				}
				else {
					std::cout << "label=<" << name << "<SUB>i</SUB>" << '>';
				}
				std::cout << ",group=\"" << name << '"';
			}
			catch (const std::out_of_range&) {
				std::cout << "label=<<I>var_cp</I>>";
			}
			std::cout << "]" << std::endl;
		}

		else {//Temporary & other											//Temporary & other
			std::cout << "[label=\"" << node->op_name_ << "\"]" << std::endl;
		}
	}

	std::string _get_cxx_node_id(const ScalarNodePtr& node) const {
		return _get_dot_node_id(node);
	}

	std::string _get_cxx_input_ids(const ScalarNodePtr& node) const {
		std::string result;
		for (uint in = 0; in < node->n_inputs_; ++in) {
			result += _get_cxx_node_id(node->inputs_[in]);
			if (in < node->n_inputs_ - 1) {
				result += ", ";
			}
		}
		return result;
	}

	std::string _get_cxx_node_definition(const ScalarNodePtr& node, const CXXVarMap& map) const {
		std::ostringstream result;

		result << "ScalarNodePtr " << _get_cxx_node_id(node) << " = ";
		switch (node->result_storage)
		{
		case ResultStorage::constant: {
			if (map.count(node->values_) == 1){
				result << "ScalarNode::create_const(node_map[\"" << map.at(node->values_) << "\"]);\n";
			}
			else {
				result << "ScalarNode::create_const(" << *node->values_ << ");\n";
			}
			break;
		}
		case ResultStorage::constant_bool: {
			if (map.count(node->values_) == 1) {
				result << "ScalarNode::create_const_bool(node_map[\"" << map.at(node->values_) << "\"]);\n";
			}
			else {
				result << "ScalarNode::create_const_bool(" << *node->values_ << ");\n";
			}
			break;
		}
		case ResultStorage::value: {
			result << "ScalarNode::create_value(node_map[\"" << map.at(node->values_) << "\"]);\n";
			break;
		}
		case ResultStorage::value_copy: {
			result << "ScalarNode::create_val_copy(node_map[\"" << map.at(node->values_) << "\"]);\n";
			break;
		}
		case ResultStorage::expr_result:
		case ResultStorage::temporary: {
			result << "ScalarNode::create<_" << node->op_name_ << "_>(" << _get_cxx_input_ids(node) << ");\n";
			break;
		}
		/*case ResultStorage::expr_result: {
			result << "ScalarNode::create_result(" << _get_cxx_node_id(node->inputs_[0]) << ", " << "???" << ");\n";
			break;
		}*/
		default:
			break;
		}
		return result.str();
	}

	std::string _get_cxx_result(const ScalarNodePtr& node, const CXXVarMap& map) const {
		return "ScalarNode::create_result(" + _get_cxx_node_id(node) + ", node_map[\"" + map.at(node->values_) + "\"]);\n";
	}

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
			ScalarNodePtr  node = nodes[i];
			for(uint in=0; in < node->n_inputs_; ++in)  {
				ScalarNodePtr  other = node->inputs_[in];
				if (other->result_idx_ != -2) {
					//BP_ASSERT(other->result_idx_ == -1);
					BP_ASSERT(other->result_storage != expr_result);
					other->result_idx_ = -2;
					nodes.push_back(other);
				}
			}
		}


		// set result_idx_ of constant nodes
		uint i_storage = 0;
		for(ScalarNodePtr  node : nodes)
			if (node->result_storage == constant || node->result_storage == constant_bool)
				node->result_idx_ = i_storage++;
		constants_end = i_storage;
		// set result_idx_ of value/result nodes
		for(ScalarNodePtr  node : nodes)
			if (node->result_storage == value || node->result_storage == expr_result)
				node->result_idx_ = i_storage++;
		values_end = i_storage;
		// set result_idx_ of value/result nodes
		for(ScalarNodePtr  node : nodes)
			if (node->result_storage == value_copy)
				node->result_idx_ = i_storage++;
		values_copy_end = i_storage;
		temp_end = i_storage; // still empty
	}


	void _topological_sort() {

		// in-degree of nodes (number of dependent nodes).
		for(ScalarNodePtr node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNodePtr node : nodes) {
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
			ScalarNodePtr  node = stack.back();
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
		for(ScalarNodePtr node : nodes) node->n_dep_nodes_ = 0;
		for(ScalarNodePtr node : nodes)
			for(uint in=0; in < node->n_inputs_; ++in) node->inputs_[in]->n_dep_nodes_ += 1;

		// Mimic expression evaluation, reversed topological order.
		for(auto it=sorted.rbegin(); it != sorted.rend(); ++it) {
			ScalarNodePtr  node = *it;
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
	void _allocate_storage(ScalarNodePtr node) {
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


	void _deallocate_storage(ScalarNodePtr node) {
		if (node->result_storage == temporary) {
			storage[node->result_idx_ - temp_end] = 0;
		}
	}










};






} // namespace details
} // namespace bparser




#endif /* INCLUDE_EXPRESSION_DAG_HH_ */
