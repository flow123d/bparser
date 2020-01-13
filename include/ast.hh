#ifndef INCLUDE_AST_HH_
#define INCLUDE_AST_HH_

#include <list>
#include <string>
#include <boost/variant.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/spirit/home/support/attributes.hpp>
#include <boost/spirit/home/qi/domain.hpp>
#include "expr.hh"



namespace bparser {

namespace ast {

typedef expr::Array (*unary_fn)(const expr::Array &);
typedef expr::Array (*binary_fn)(const expr::Array &, const expr::Array &);


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
 * Phoenix Helpers
 */

// unary expression factory
struct make_unary_f {
	unary_op operator()(unary_fn op, operand const& lhs) const {
		return {op, lhs};
	}
};
boost::phoenix::function<make_unary_f> make_unary;


// binary expression factory
struct make_binary_f {
	binary_op operator()(binary_fn op, operand const& lhs, operand const& rhs) const {
		return {op, lhs, rhs};
	}
};
::boost::phoenix::function<make_binary_f> make_binary;

// assign expression factory
struct make_assign_f {
	assign_op operator()(std::string lhs, operand const& rhs) const {
		return {lhs, rhs};
	}
};
::boost::phoenix::function<make_assign_f> make_assign;


expr::Array semicol_fn(const expr::Array &a, const expr::Array &b) {
	return b;
}

} // namespace ast

} // namespace bparser







#endif //INCLUDE_PARSER_HH_
