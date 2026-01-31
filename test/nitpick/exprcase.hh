/*
 * exprcase.hh
 * Holds the memory and other config for generating and running DAGs
 *
 *  Created on: Jan 20, 2026
 *      Author: LV
 */


#ifndef EXPRCASE_HH
#define EXPRCASE_HH

#include "arena_resource.hh"
#include "parser.hh"

using namespace bparser;

namespace bparser {

	typedef std::unordered_map<std::string, double*> NodeMap;

	class ExprCase {
	public:
		enum VarType {
			Variable,
			VarCopy,
			Const
		};

		using string = std::string;
		using Shape = bparser::Shape;
		typedef std::pair<VarType, Shape> VarInfo;



	public: //Constructors
		ExprCase(uint max_vec_size, void* buffer, size_t buffer_size, size_t simd_size) :
			max_vec_size(max_vec_size),
			arena(std::make_shared<PatchArena>(buffer, buffer_size, simd_size)),
			parser(Parser(max_vec_size))
		{
			;
		}

		ExprCase(uint max_vec_size, void* buffer, size_t buffer_size) : 
			ExprCase(max_vec_size, buffer, buffer_size, bparser::get_simd_size())
		{ 
			;
		}
	protected:

		uint max_vec_size;

		PatchArenaPtr arena;
		Parser parser;
		std::string parseExpr;

		std::unordered_map<string, VarInfo> variables;
		std::unordered_map<string, std::vector<double>> variable_consts;

		std::unordered_map<double*, std::string> inv_map;  //given to the print_in_cxx function to generate the file
		NodeMap node_map; //is used by the generated file to use the same double pointers as the parser
		

	public: //Methods for the *_def.hh files
		void set_variable(string name, Shape shape){
			variables[name] = VarInfo(Variable, shape);
		}
		void set_var_copy(string name, Shape shape){
			variables[name] = VarInfo(VarCopy, shape);
		}
		void set_constant(string name, Shape shape, std::vector<double> values){
			variables[name] = VarInfo(Const, shape);
			variable_consts[name] = values;
		}
		void set_result_shape(Shape shape) {
			variables["_result_"] = VarInfo(Variable, shape);
		}
		void parse(string expression) {
			parseExpr = expression;

			parser.parse(expression);
		}

	public: //Methods for the gen/run files

		void deallocate() {
			arena->reset();
		}

		void allocate(double vec_size) {
			for (const auto& [varName, varInfo] : variables) {
				VarType t = varInfo.first;
				Shape shape = varInfo.second;
				double n_values = numel(shape);
				double* ptr = arena->allocate_simd<double>(n_values * vec_size);

				if (t == Const) {
					map_const(varName, variable_consts[varName], shape);
				}
				else {
					map_variable(varName, ptr, shape);
				}
				

				if (t == Variable) {
					parser.set_variable(varName, shape, ptr);
				}else if (t == VarCopy) {
					parser.set_var_copy(varName, shape, ptr);
				}
				else if (t == Const) {
					parser.set_constant(varName, shape, variable_consts[varName]);
				}

			}
		}

		PatchArenaPtr get_patch_arena() const {
			return arena;
		}
		string get_expression() const {
			return parseExpr;
		}

		const std::unordered_map<double*, std::string>& get_inv_map() const{
			return inv_map;
		}
		const NodeMap& get_node_map() const{
			return node_map;
		}

		Parser& get_parser() {
			return parser;
		}

	protected: //Helper methods

		//Put the variable names and ptrs in the maps for gen/run use
		void map_variable(string name,double* pointer, const Shape& shape) {
			Array array = Array::value(pointer, max_vec_size, shape);

			for (MultiIdx idx(array.range()); idx.valid(); idx.inc_src()) {

				std::string var_name = get_var_name(name, idx.indices());
				inv_map[array[idx]->values_] = var_name;
				node_map[var_name] = array[idx]->values_;
			}
		}

		void map_const(string name, const std::vector<double>& values, const Shape& shape) {
			Array array = Array::constant(values, shape); //useless

			for (MultiIdx idx(array.range()); idx.valid(); idx.inc_src()) {

				std::string var_name = get_var_name(name, idx.indices());
				//inv_map[array[idx]->values_] = var_name;
				//node_map[var_name] = array[idx]->values_;
			}
		}

		inline double numel(const Shape& s) {
			double res = 1;
			for (const uint& n : s) {
				res *= n;
			}
			return res;
		}

		// v1, (1,2,3) => "v1[1,2,3]"
		std::string get_var_name(const std::string& name, const bparser::MultiIdx::VecUint& indices) {
			bparser::MultiIdx::VecUint::size_type size(indices.size());

			std::ostringstream result;
			result << name;
			result << "[";
			for (bparser::MultiIdx::VecUint::size_type i = 0; i < size; i++) {
				result << indices.at(i);
				if (i != size - 1) {
					result << ',';
				}
			}
			result << "]";
			return result.str();
		}
	};
}

#endif //EXPRCASE_HH