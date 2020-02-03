#ifndef INCLUDE_AST_HH_
#define INCLUDE_AST_HH_

#include <list>
#include <string>
#include <boost/variant.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/spirit/home/qi/domain.hpp>

#include "array.hh"



namespace bparser {

namespace ast {
/**
 * Abstract Syntax Tree resulting from the parsed grammar.
 * - Various helper structs and classes for AST construction (via. BOOST functional programming tools)
 * - Boost Variant type visitors for AST processing.
 * - Finally the AST is translated into DAG (direct acyclic graph) of ScalarNode operations
 *   grouped into vector structures via. Array classes.
 */




struct unary_fn {
	typedef Array (*type)(const Array &);
	std::string repr;
	type fn;
};

struct binary_fn {
	typedef Array (*type)(const Array &, const Array &);
	std::string repr;
	type fn;
};




struct nil {};
struct unary_op;
struct binary_op;
struct assign_op;

// clang-format off
typedef boost::variant<
        nil // can't happen!
        , double
        , std::string
        , boost::recursive_wrapper<unary_op>
        , boost::recursive_wrapper<binary_op>
		, boost::recursive_wrapper<assign_op>
        >
operand;
// clang-format on

struct unary_op {
    unary_fn op;
    operand rhs;
};

struct binary_op {
    binary_fn op;
    operand lhs;
    operand rhs;
};


struct assign_op {
    std::string lhs;
    operand rhs;
};

/**
 * operand visitors
 */


struct print_vis : public boost::static_visitor<> {
	/**
	 * Print AST.
	 */
	mutable std::stringstream ss;

    explicit print_vis()
    {}


    void operator()(nil x) const {
    	ss << "NULL";
    }

    void operator()(double x) const
    {ss << x; }

    void operator()(std::string const &x) const  {
    	ss << "<" << x << ">";
    }


    void operator()(unary_op const &x) const {
    	ss << x.op.repr << "(";
    	boost::apply_visitor(*this, x.rhs);
    	ss << ")";
    }

    void operator()(binary_op const &x) const {
    	ss << x.op.repr << "(";
    	boost::apply_visitor(*this, x.lhs);
    	ss << ", ";
    	boost::apply_visitor(*this, x.rhs);
    	ss << ")";
    }

    void operator()(assign_op const &x) const  {
    	ss << x.lhs << " = ";
    	boost::apply_visitor(*this, x.rhs);
    }


};

/**
 * Print AST function
 */
inline std::string print(operand const& ast_root) {
	print_vis pv;
	boost::apply_visitor(pv, ast_root);
	return pv.ss.str();
}

/**
 * Lazy AST print.
 */
struct print_ast_t {
	typedef bool result_type;

	bool operator()(operand const& x) const {
		print(x);
		return true;
	}
};
inline boost::phoenix::function<print_ast_t> lazy_print;



// unary expression factory
struct make_unary_f {
	unary_op operator()(unary_fn op, operand const& lhs) const {
		// std::cout << "make_unary: " << lhs.which() << "\n";
		return {op, lhs};
	}
};
inline boost::phoenix::function<make_unary_f> make_unary;


// binary expression factory
struct make_binary_f {
	binary_op operator()(binary_fn op, operand const& lhs, operand const& rhs) const {
		// std::cout << "make_binary: " << lhs.which() << ", " << rhs.which() << "\n";
		return {op, lhs, rhs};
	}
};
inline boost::phoenix::function<make_binary_f> make_binary;

struct make_relational_f {
	binary_op operator()(binary_fn op, operand const& lhs, operand const& rhs, operand const& chained) const {
		std::cout << "make_binary: " << lhs.which() << ", " << rhs.which() << ", " << print(chained) << "\n";
		return {op, lhs, rhs};
	}
};
inline boost::phoenix::function<make_relational_f> make_relational;



// assign expression factory
struct make_assign_f {
	assign_op operator()(std::string lhs, operand const& rhs) const {
		return {lhs, rhs};
	}
};
inline boost::phoenix::function<make_assign_f> make_assign;


inline Array semicol_fn(const Array &a, const Array &b) {
	return b;
}




/**
 * Construct ScalarNode operation DAG via. Arrays to expand vector operations
 * into scalar operations.
 */
struct make_array {
    typedef Array result_type;
    mutable std::map<std::string, Array> symbols;

    explicit make_array(std::map<std::string, Array> const &symbols)
    : symbols(symbols)
    {}

    result_type operator()(nil) const {
        BOOST_ASSERT(0);
        return Array();
    }

    result_type operator()(double x) const
    {
    	return Array::constant({x});
    }

    result_type operator()(std::string const &x) const  {
        auto it = symbols.find(x);
        if (it != symbols.end()) {
        	return it->second;
        } else {
        	std::ostringstream s;
        	s << "Undefined var: " << x << "\n";
        	// We do not call visitor for the assign_op's 'lhs' so this must be error.
        	Throw(s.str());
        }
    }

    result_type operator()(unary_op const &x) const {
        return x.op.fn(boost::apply_visitor(*this, x.rhs));
    }

    result_type operator()(binary_op const &x) const {
    	result_type lhs = boost::apply_visitor(*this, x.lhs);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
        return x.op.fn(lhs, rhs);
    }

    result_type operator()(assign_op x) const  {
    	//std::string& var_name = boos);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	symbols[x.lhs] = rhs;
        return rhs;
    }

};




/**
 * Extract free symbols/variables.
 */
struct get_variables {
    typedef std::vector<std::string> result_type;

    static result_type merge(result_type a, result_type b) {
    	result_type res(a.begin(), a.end());
    	res.insert(res.end(), b.begin(), b.end());
    	return res;
    }

    explicit get_variables()
    {}


    result_type operator()(nil) const {
        BOOST_ASSERT(0);
        return {};
    }

    result_type operator()(double x) const
    { return {}; }

    result_type operator()(std::string const &x) const  {
        auto it = expr_defs.find(x);
        if (it == expr_defs.end()) {
        	return {x};
        }
        return {};
    }


    result_type operator()(unary_op const &x) const {
        return boost::apply_visitor(*this, x.rhs);
    }

    result_type operator()(binary_op const &x) const {
    	result_type lhs = boost::apply_visitor(*this, x.lhs);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
        return merge(lhs, rhs);
    }

    result_type operator()(assign_op const &x) const  {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	expr_defs.insert(x.lhs);
        return rhs;
    }



private:
    mutable std::set<std::string> expr_defs;
};





/**
 * Remove 'nil' nodes from AST. Should not be necessary.
 */
struct remove_nil {
    typedef operand result_type;


    explicit remove_nil()
    {}


    result_type operator()(nil) const {
        return nil();
    }

    result_type operator()(double x) const
    { return x; }

    result_type operator()(std::string const &x) const  {
        return x;
    }


    result_type operator()(unary_op const &x) const {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	if (rhs.type() != typeid(nil)) return x;
    	return rhs;
    }

    result_type operator()(binary_op const &x) const {
    	result_type lhs = boost::apply_visitor(*this, x.lhs);
    	if (lhs.type() != typeid(nil)) return x;
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	if (rhs.type() != typeid(nil)) return x;
    	return rhs;
    }

    result_type operator()(assign_op const &x) const  {
    	BP_ASSERT(x.lhs.size() > 0);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	BP_ASSERT(rhs.type() != typeid(nil));
        return x;
    }



};


} // namespace ast

} // namespace bparser






#endif //INCLUDE_PARSER_HH_
