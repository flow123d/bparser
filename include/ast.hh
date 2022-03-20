#ifndef INCLUDE_AST_HH_
#define INCLUDE_AST_HH_

#include <list>
#include <string>
#include <algorithm>

#include <boost/variant.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/spirit/home/qi/domain.hpp>

#include "config.hh"
#include "array_ast_interface.hh"

/**
 * TODO:
 * 1. We can possibly skip the AST step and directly create ScalareNode DAG through the Array interface.
 *    This needs some modifications:
 *    a) We have to find a way to have PHOENIX delayed functions replaced by delayed methdos of some object.
 *       We need one class (PrintBackend) to print AST of the grammar.
 *       We need another class (EvalBackend) to construct the evaluation DAG through the Arrays.
 *    b) Delayed construction of ValueNodes from identifiers. The `get_variables` step has to be replaced
 *       by storing uncomplete ValueNodes as a map in EvalBackend.
 *
 *    Advantage: less code, possible small speed gain in parsing phase
 *    Disadvantage: no mean for AST processing, however AST is a bit arbitrary structure anyway
 */

namespace bparser {

namespace ast {
/**
 * Abstract Syntax Tree resulting from the parsed grammar.
 * - Various helper structs and classes for AST construction (via. BOOST functional programming tools)
 * - Boost Variant type visitors for AST processing.
 * - Finally the AST is translated into DAG (direct acyclic graph) of ScalarNode operations
 *   grouped into vector structures via. Array classes.
 */


//struct unary_fn {
//	typedef Array (*function_type)(const Array &);
//	std::string repr;
//	function_type fn;
//};
//
//struct binary_fn {
//	typedef Array (*function_type)(const Array &, const Array &);
//	std::string repr;
//	function_type fn;
//};
//
//struct ternary_fn {
//	typedef Array (*function_type)(const Array &, const Array &, const Array &);
//	std::string repr;
//	function_type fn;
//};



//struct nil {};
struct list;
struct call;
struct assign_op;

// clang-format off
typedef boost::variant<
        //nil, // can't happen!  TODO: remove
        double
		, int
        , std::string
        , boost::recursive_wrapper<list>
        , boost::recursive_wrapper<call>
		, boost::recursive_wrapper<assign_op>
        >
operand;
// clang-format on

/// a function: Array -> Array
struct list {
    operand head; // list; any non-list type marks end of recursion; none-int is convention
    operand item; //
};


/// a function: (Array, Array) -> Array
struct call {
    NamedArrayFn op;
    operand arg_list;
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

    bool print_list(operand x) const {
    	if (list *l = boost::get<list>(&x)) {
			if (print_list(l->head))
				ss << ",";
			boost::apply_visitor(*this, l->item);
			return true;
    	}
    	return false;
    }

    void operator()(double x) const
    {ss << x; }

    void operator()(int x) const
    {
    	if (x == none_int) {
    		ss << "none_idx";
    	} else {
    		ss << x;
    	}
    }

    void operator()(std::string const &x) const  {
    	ss << "`" << x << "`";
    }

    void operator()(call const &x) const {
    	ss << x.op.repr;
    	ss << '(';
    	//std::cout << x.op.repr << "\n";
    	//assert_list(x.arg_list);
    	print_list(x.arg_list);
    	ss << ')';

    }

    void operator()(list x) const {
    	ss << '[';
    	print_list(x);
    	ss << ']';
    	//BP_ASSERT(false);
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

inline std::ostream &operator<<(std::ostream &os, const operand& a) {
	print_vis pv;
	boost::apply_visitor(pv, a);
    os << pv.ss.str();
    return os;
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
// this and following declarations can be used with C++17
//inline boost::phoenix::function<print_ast_t> lazy_print;
BOOST_PHOENIX_ADAPT_CALLABLE(lazy_print, print_ast_t, 1)






/**
 * Construct ScalarNode operation DAG via. Arrays to expand vector operations
 * into scalar operations.
 */
struct make_array {
    typedef ParserResult result_type;

    mutable std::map<std::string, Array> symbols;

    explicit make_array(std::map<std::string, Array> const &symbols)
    : symbols(symbols)
    {}


    ResultList make_list(operand x) const {
    	ResultList alist;
    	if (list *l = boost::get<list>(&x)) {
    		alist = make_list(l->head);
			ParserResult aitem = boost::apply_visitor(*this, l->item);
			alist.push_back(aitem);
    	}
    	return alist;
    }


//    result_type operator()(nil) const {
//        BOOST_ASSERT(0);
//        return Array();
//    }

    result_type operator()(double x) const
    {
    	return x;
    }

    result_type operator()(int x) const
    {return x;}


    result_type operator()(std::string const &x) const  {
        auto it = symbols.find(x);
        BP_ASSERT(it != symbols.end());
        if (! it->second.is_none()) {
        	return it->second;
        } else {
        	// We do not call visitor for the assign_op's 'lhs' so this must be error.
        	Throw() << "Undefined var: " << x << "\n";
        }
    }

    result_type operator()(call const &x) const {
    	//assert_list(x.arg_list);
    	ResultList al = make_list(x.arg_list);
    	//std::cout << "apply fn:" << x.op.repr << "\n";
    	result_type res = boost::apply_visitor(call_visitor(al), x.op.fn);
    	return res;
    }

    result_type operator()(list x) const {
    	ResultList al = make_list(x);
    	if (std::all_of(al.begin(), al.end(), [](result_type x){return check_type<int>(x);})) {
        	ListInt alist;
        	for(ParserResult item : al)
        		alist.push_back(get_type<int>(item));
    		return alist;
    	}

    	std::vector<Array> alist;
    	for(ParserResult item : al)
    		alist.push_back(get_array(item));
    	return Array::stack_zero(alist);
    }


    result_type operator()(assign_op x) const  {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	if (Array* rhs_array = boost::get<Array>(&rhs)) {
    		symbols[x.lhs] = *rhs_array;
    		return rhs;
    	}
    	Throw() << "Internal error.\n"
    		<< "Wrong type: " << "operand:" << rhs.which() << "\n"
			<< "Expected: Array\n";
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

    result_type merge_list(operand x) const {
    	if (list *l = boost::get<list>(&x)) {
        	result_type head_vars = merge_list(l->head);
        	result_type item_vars = boost::apply_visitor(*this, l->item);
        	return merge(head_vars, item_vars);
    	} else {
    		return result_type();
    	}
    }


//    result_type operator()(nil) const {
//        BOOST_ASSERT(0);
//        return {};
//    }

    result_type operator()(double UNUSED(x)) const
    { return {}; }

    result_type operator()(std::string const &x) const  {
        auto it = expr_defs.find(x);
        if (it == expr_defs.end()) {
        	return {x};
        }
        return {};
    }


    result_type operator()(call const &x) const {
    	//assert_list(x.arg_list);
    	return merge_list(x.arg_list);
    }

    result_type operator()(list x) const {
    	return merge_list(x);
    }

    result_type operator()(assign_op const &x) const  {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	expr_defs.insert(x.lhs);
        return rhs;
    }



private:
    mutable std::set<std::string> expr_defs;
};





// Grammar PHOENIX lazy functions
// using BOOST_PHOENIX_ADAPT_CALLABLE (instead of ..._ADAPT_FUNCTION)
// in order to avoid specification of the return type.

// unary expression factory
struct make_call_f {
	call operator()(NamedArrayFn op, operand const& first) const {
		return {op, first};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_call, make_call_f, 2)


// unary expression factory
struct make_unary_f {
	call operator()(NamedArrayFn op, operand const& first) const {
		// std::cout << "make_unary: " << first.which() << "\n";
		list l = {0.0, first};
		return {op, l};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_unary, make_unary_f, 2)

// nulary expression factory
struct make_const_f {
	// In order to minimize number of differrent operations in AST
	// the nulary 'fn' must in fact accept single double argument.
	call operator()(std::string repr, ArrayFnUnary fn) const {
		NamedArrayFn nulary_fn = {repr, fn};
		list l = {0.0, 0.0};
		return {nulary_fn, l};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_const, make_const_f, 2)


// binary expression factory
struct make_binary_f {
	call operator()(NamedArrayFn op, operand const& first, operand const& second) const {
		// std::cout << "make_binary: " << first.which() << ", " << second.which() << "\n";
		list l1 = {0.0, first};
		list l2 = {l1, second};
		return {op, l2};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_binary, make_binary_f, 3)

struct make_ternary_f {
	call operator()(NamedArrayFn op, operand const& first, operand const& second, operand const& third) const {
		// std::cout << "make_binary: " << first.which() << ", " << second.which() << "\n";
		list l1 = {0.0, first};
		list l2 = {l1, second};
		list l3 = {l2, third};
		return {op, l3};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_ternary, make_ternary_f, 4)



struct make_chained_comparison_f {
	call operator()(NamedArrayFn op, operand const& head, operand const& other) const {
		//std::cout << "make_binary: " << first.which() << ", " << second.which() << ", " << print(chained) << "\n";
		list l1 = {0.0, head};
		list l2 = {l1, other};
		return {op, l2};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_chained_comparison, make_chained_comparison_f, 3)


struct make_list_f {
	operand operator()(operand const& head, operand const& item) const {
		return list({head, item});
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_list, make_list_f, 2)


// assign expression factory
struct make_assign_f {
	assign_op operator()(std::string lhs, operand const& rhs) const {
		return {lhs, rhs};
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(make_assign, make_assign_f, 2)

// Convert 'boost::optional' to 'operand'.
struct treat_optional_f {
	operand operator()(boost::optional<operand> const& v, operand const& default_) const {
		if (v == boost::none) {
			return default_;
		} else {
			return *v; // extract optional value
		}
	}
};
BOOST_PHOENIX_ADAPT_CALLABLE(treat_optional, treat_optional_f, 2)



inline Array semicol_fn(const Array & UNUSED(a), const Array &b) {
	return b;
}








} // namespace ast

} // namespace bparser






#endif //INCLUDE_PARSER_HH_
