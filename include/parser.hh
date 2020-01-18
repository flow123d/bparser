/*
 * parser.hh
 *
 *  Created on: Jan 11, 2020
 *      Author: jb
 */


#ifndef INCLUDE_PARSER_HH_
#define INCLUDE_PARSER_HH_


#include "ast.hh"
#include "evaluator.hh"
#include "processor.hh"
#include "grammar.hh"
#include "expr.hh"


#include <map>
#include <memory>
#include <string>

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
	std::map<std::string, expr::Array> symbols_;
	std::vector<std::string> free_variables;
	expr::Array result_array;
	Processor * processor;

public:
    /// @brief Constructor
    Parser(uint max_vec_size)
	: max_vec_size(max_vec_size), processor(nullptr)
	{}

    /// @brief Destructor
    ~Parser() {}

    void destroy_processor() {
    	if (processor != nullptr) processor->~Processor();
    }
    /// @brief Parse the mathematical expression into an abstract syntax tree
    ///
    /// @param[in] expr The expression given as a std::string
    void parse(std::string const &expr) {
        ast::operand ast_;

        std::string::const_iterator first = expr.begin();
        std::string::const_iterator last = expr.end();

        boost::spirit::ascii::space_type space;
        bool r = qi::phrase_parse(
            first, last, grammar(), space,
            ast_);

        if (!r || first != last) {
            std::string rest(first, last);
            Throw("Parsing failed at " + rest); // NOLINT
        }
        //std::cout << "Parsing OK. : " << "\n";
        //print(ast_);

        //ast = boost::apply_visitor(ast::remove_nil(), ast_);
        ast = ast_;
        //ASSERT(ast.type() != typeid(ast::nil));
        free_variables  = boost::apply_visitor(ast::get_variables(), ast);
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
    std::vector<std::string> const &variables() {
    	return free_variables;
    }

    /**
     * Set given name to be a variable of given shape with values at
     * given address 'variable_space'.
     *
     */
    void set_variable(std::string name, std::vector<uint> shape, double *variable_space) {
    	symbols_[name] = expr::Array::value(variable_space, max_vec_size, shape);
    }

    /**
     * Set given name to be a constant of given shape with flatten values
     * given by the 'const_value' vector.
     *
     */
    void set_constant(std::string name, std::vector<uint> shape, std::vector<double> const_value) {
    	symbols_[name] = expr::Array::constant(const_value, shape);
    }

    /// @brief Create processor of the expression from the AST.
    ///
    /// All variable names have to be set before this call.
    /// TODO: set result variable
    void compile() {
    	auto res_it = symbols_.find("_result_");
    	if (res_it == symbols_.end())
    		Throw("No '_result_' set.");

        expr::Array res_array = boost::apply_visitor(ast::make_array(symbols_), ast);
        result_array = res_array.make_result(res_it->second);

        destroy_processor();
		ScalarExpression se(result_array.elements());

		se.print_in_dot();
		processor = Processor::create_processor_(se, max_vec_size);
    }


    /// @brief Set new subset of the 'max_vec_size' vectors.
    /// Only this subset is evuluated by the processor.
    void set_subset(std::vector<uint> const &subset) {
    	ASSERT(processor != nullptr);
    	processor->set_subset(subset);
    }

    void run() {
    	processor->run();
    }


};


} // namespace bparser




#endif /* INCLUDE_PARSER_HH_ */
