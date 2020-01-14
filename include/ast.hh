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

struct unary_fn {
	typedef expr::Array (*type)(const expr::Array &);
	std::string repr;
	type fn;
};

struct binary_fn {
	typedef expr::Array (*type)(const expr::Array &, const expr::Array &);
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
 * Phoenix Helpers
 */





struct print_vis : public boost::static_visitor<> {
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


std::string print(operand const& x) {
	print_vis pv;
	boost::apply_visitor(pv, x);
	return pv.ss.str();
}

struct print_ast_t {
	typedef bool result_type;

	bool operator()(operand const& x) const {
		print(x);
		return true;
	}
};

boost::phoenix::function<print_ast_t> lazy_print;



// unary expression factory
struct make_unary_f {
	unary_op operator()(unary_fn op, operand const& lhs) const {
		// std::cout << "make_unary: " << lhs.which() << "\n";
		return {op, lhs};
	}
};
boost::phoenix::function<make_unary_f> make_unary;


// binary expression factory
struct make_binary_f {
	binary_op operator()(binary_fn op, operand const& lhs, operand const& rhs) const {
		// std::cout << "make_binary: " << lhs.which() << ", " << rhs.which() << "\n";
		return {op, lhs, rhs};
	}
};
::boost::phoenix::function<make_binary_f> make_binary;

struct make_relational_f {
	binary_op operator()(binary_fn op, operand const& lhs, operand const& rhs, operand const& chained) const {
		std::cout << "make_binary: " << lhs.which() << ", " << rhs.which() << ", " << print(chained) << "\n";
		return {op, lhs, rhs};
	}
};
::boost::phoenix::function<make_relational_f> make_relational;



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
