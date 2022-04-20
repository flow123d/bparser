/*
 * parser.hh
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */


#ifndef INCLUDE_PARSER_HH_
#define INCLUDE_PARSER_HH_


#include <map>
#include <memory>
#include <string>
#include <vector>
#include "array.hh"
#include "ast.hh"
#include "processor.hh"
#include "grammar.hh"

namespace bparser {

/// @brief Parse a mathematical expression
///
/// This can parse and evaluate a mathematical expression for a given
/// symbol table using Boost.Spirit X3.  The templates of Boost.Spirit
/// are very expensive to parse and instantiate, which is why we hide
/// it behind an opaque pointer.
///
/// The drawback of this approach is that calls can no longer be
/// inlined and because the pointer crosses translation unit
/// boundaries, dereferencing it can also not be optimized out at
/// compile time.  We have to rely entirely on link-time optimization
/// which might be not as good.
///
/// The pointer to the implementation is a std::unique_ptr which makes
/// the class not copyable but only moveable.  Copying shouldn't be
/// required but is easy to implement.



class Parser {

	ast::operand ast;
	uint max_vec_size;
    uint simd_size;
	std::map<std::string, Array> symbols_;
	Array result_array_;
	ProcessorBase * processor;
	std::vector<double> tmp_result;

public:
    /** @brief Constructor
     * max_vec_size - size of single array component in doubles
     */
    Parser(uint max_vec_size, uint simd_size)
	: max_vec_size(max_vec_size), simd_size(simd_size), processor(nullptr), tmp_result()
	{}
    // Parser(uint max_vec_size, uint simd_size)
	// : max_vec_size(max_vec_size), simd_size(simd_size), processor(nullptr), tmp_result(nullptr)
	// {}

    /// @brief Destructor
    ~Parser() {
    	destroy_processor();
    }

    void destroy_processor() {
        // arena_.destroy();
    	// if (tmp_result != nullptr) delete [] tmp_result;
    	if (processor != nullptr) processor->~ProcessorBase();
    	processor = nullptr;
    }
    /// @brief Parse the mathematical expression into an abstract syntax tree
    ///
    /// @param[in] expr The expression given as a std::string
    void parse(std::string const &expr) {
    	parse_expr(expr, ast);

        //std::cout << "Parsing OK. : " << "\n";
        //std::cout << ast::print(ast) << "\n";

        //ASSERT(ast.type() != typeid(ast::nil));

    	std::vector<std::string> free_variables;
    	free_variables = boost::apply_visitor(ast::get_variables(), ast);
    	for(std::string &s: free_variables) {
        	symbols_[s] = Array(); // none array

        }


        // default constants
        set_constant("e", {}, {boost::math::constants::e<double>()});
        set_constant("pi", {}, {boost::math::constants::pi<double>()});
        //set_constant("phi", {}, {boost::math::constants::phi<double>()});
        // rounding precision
        set_constant("epsilon", {}, {std::numeric_limits<double>::epsilon()});

        //_optimize();
    }

    std::string print_ast() {
    	return ast::print(ast);
    }



    void _optimize() {
    	//ast = boost::apply_visitor(ast::ConstantFolder(), ast);
    }

    /**
     * @brief Return names of known symbols.
     */
    std::vector<std::string> symbols() const {
    	std::vector<std::string> keys;
    	for(auto s : symbols_)
    		keys.push_back(s.first);
    	return keys;
    }

    /**
     * @brief Return names (undefined) variables in the expression.
     */
    std::vector<std::string> free_symbols() {
    	std::vector<std::string> keys;
    	for(auto s : symbols_) {
    		if (s.second.is_none())	keys.push_back(s.first);
    	}
    	return keys;
    }

    /**
     * Set given name to be a variable of given shape with values at
     * given address 'variable_space'.
     *
     * Unused variables and constants are ignored.
     *
     */
    void set_variable(std::string name, std::vector<uint> shape, double *variable_space) {
    	symbols_[name] = Array::value(variable_space, max_vec_size, shape);
    }

    /**
     * Set given name to be a constant of given shape with flatten values
     * given by the 'const_value' vector.
     *
     */
    void set_constant(std::string name, std::vector<uint> shape, std::vector<double> const_value) {
    	symbols_[name] = Array::constant(const_value, shape);
    }

    /// @brief Create processor of the expression from the AST.
    ///
    /// All variable names have to be set before this call.
    /// TODO: set result variable
    // ExpressionDAG compile() {
    void compile() {
    	destroy_processor();

        ParserResult res_array = boost::apply_visitor(ast::make_array(symbols_), ast);

        Array array = get_array(res_array);
		Shape result_shape = array.shape();
		auto res_it = symbols_.find("_result_");
		if (res_it == symbols_.end()) {
			// TODO: replace by storing result in the temporary variable of the processor
			// tmp_result = new double[shape_size(result_shape) * max_vec_size];
            tmp_result.resize(shape_size(result_shape) * max_vec_size);
			result_array_ = Array::value(&tmp_result[0], max_vec_size, result_shape);
			result_array_ = array.make_result(result_array_);
		} else {
			result_array_ = array.make_result(res_it->second);
		}

		ExpressionDAG se(result_array_.elements());
        // return se;
		//se.print_in_dot();
		processor = create_processor(se, max_vec_size, simd_size);
    }

    Array result_array() {
    	return result_array_;
    }

    double * tmp_result_ptr() {
    	return &tmp_result[0];
    }

    /// @brief Set new subset of the 'max_vec_size' vectors.
    /// Only this subset is evuluated by the processor.
    void set_subset(std::vector<uint> const &subset) {
    	BP_ASSERT(processor != nullptr);
    	processor->set_subset(subset);
    }

    void run() {
    	processor->run();
    }

    void set_processor(ProcessorBase *p)
    {
        processor = p;
    }


};


} // namespace bparser



#endif /* INCLUDE_PARSER_HH_ */
