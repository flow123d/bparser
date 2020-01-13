#ifndef INCLUDE_GRAMMAR_HH_
#define INCLUDE_GRAMMAR_HH_


#include <iostream>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/phoenix.hpp>

//#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#include <boost/spirit/include/qi.hpp>


#include "ast.hh"
#include "scalar_expr.hh"
#include "expr.hh"
//#include "math.hpp"



namespace bparser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


namespace parser {




template<class T>
ast::unary_fn unary_array() {
	return static_cast<ast::unary_fn>(&(expr::Array::unary_op<T>));
}

template<class T>
ast::binary_fn binary_array() {
	return static_cast<ast::binary_fn>(&(expr::Array::binary_op<T>));
}



struct expectation_handler {
    template <typename>
    struct result {
        typedef void type;
    };

    template <typename Iterator>
    void operator()(Iterator first, Iterator last,
                    boost::spirit::info const &info) const {
        std::stringstream msg;
        msg << "Expected " << info << " at \"" << std::string(first, last)
            << "\"";

        throw std::runtime_error(msg.str()); // NOLINT
    }
};


template <typename Iterator>
struct grammar : qi::grammar<Iterator, ast::operand(), ascii::space_type> {
    expectation_handler err_handler;


    qi::rule<Iterator, ast::operand(), ascii::space_type>
    	primary, logical, equality, relational, additive, multiplicative, factor, unary, binary, program, definition, assignment;
    qi::rule<Iterator, std::string()> variable;


    qi::symbols<typename std::iterator_traits<Iterator>::value_type, double>
        constant;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                ast::unary_fn>
        ufunc, unary_op;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                ast::binary_fn>
        bfunc, additive_op, multiplicative_op, logical_op, relational_op,
		equality_op, power, semicol_op;


    grammar()
    : grammar::base_type(program)
    {
    	using namespace qi;

//        qi::alnum_type alnum;
//        qi::alpha_type alpha;
//        qi::double_type double_;
//        qi::lexeme_type lexeme;
//        qi::raw_type raw;

        // clang-format off

        constant.add
            ("e"      , boost::math::constants::e<double>())
            ("epsilon", std::numeric_limits<double>::epsilon())
            ("phi"    , boost::math::constants::phi<double>())
            ("pi"     , boost::math::constants::pi<double>())
            ;

        ufunc.add
            ("abs"  , unary_array<_abs_>())
            ("acos" , unary_array<_acos_>())
            ("asin" , unary_array<_asin_>())
            ("atan" , unary_array<_atan_>())
            ("ceil" , unary_array<_ceil_>())
            ("cos"  , unary_array<_cos_>())
            ("cosh" , unary_array<_cosh_>())
            ("deg"  , &expr::deg_fn)
            ("exp"  , unary_array<_exp_>())
            ("floor", unary_array<_floor_>())
            ("isinf", unary_array<_isinf_>())
            ("isnan", unary_array<_isnan_>())
            ("log"  , unary_array<_log_>())
            ("log10", unary_array<_log10_>())
            ("rad"  , &expr::rad_fn)
            ("sgn"  , unary_array<_sgn_>())
            ("sin"  , unary_array<_sin_>())
            ("sinh" , unary_array<_sinh_>())
            ("sqrt" , unary_array<_sqrt_>())
            ("tan"  , unary_array<_tan_>())
            ("tanh" , unary_array<_tanh_>())
            ;

        bfunc.add
            ("atan2", binary_array<_atan2_>())
            ("pow"  , binary_array<_pow_>())
            ;

        unary_op.add
            ("+", &expr::unary_plus)
            ("-", unary_array<_minus_>())
            ("!", unary_array<_neg_>())
            ;

        additive_op.add
            ("+", binary_array<_add_>())
            ("-", binary_array<_sub_>())
            ;

        multiplicative_op.add
            ("*", binary_array<_mul_>())
            ("/", binary_array<_div_>())
            ("%", binary_array<_mod_>())
            ;

        logical_op.add
            ("&&", binary_array<_and_>())
            ("||", binary_array<_or_>())
            ;

        relational_op.add
            ("<" , binary_array<_lt_>())
            ("<=", binary_array<_le_>())
            (">" , &expr::gt_op)
            (">=", &expr::ge_op)
            ;

        equality_op.add
            ("==", binary_array<_eq_>())
            ("!=", binary_array<_ne_>())
            ;

        power.add
            ("**", binary_array<_pow_>())
            ;

//        semicol_op.add
//            (";", &ast::semicol_fn)  // TODO: possibly some special function
//			;

        /**
         * Syntax notes:
         * lhs >> rhs : match if rhs follows lhs
         * lhs > rhs  : same, but throw if there is lhs and no rhs
         * -(expr)     : optional rule
         * *(expr)     : match zero or more repetitions
         * |           : alternative (tried sequantionaly)
         *
         * Every rule have associated type it produce in  AST given as
         * second template parameter to qi::rule.
         *
         * boost::phoenix try to treat given type as the default type
         * of the associated storage (called attribute)
         *
         * For the list of rule operators and their associated attributes see:
         * https://www.boost.org/doc/libs/1_72_0/libs/spirit/doc/html/spirit/qi/quick_reference/qi_parsers/operator.html
         * For attributes introduction:
         * http://boost-spirit.com/home/articles/attribute_handling/the-magical-power-of-attributes-in-spirit-primitives/
         * ... and sequels.
         */

        /**
         * TODO: need to construct binary_op attribute, but
         * the order of opeartions doesn't match
         */
        program = -(definition > ';')[_val = _1] >> logical[_val = ast::make_binary(&ast::semicol_fn, _val,  _1)]
			;

        definition =
        		assignment[_val = _1] >> *(semicol_op > definition)[_val = ast::make_binary(_1, _val, _2)]
				;

        assignment =
        	(variable > '=' > logical)[_val = ast::make_assign(_1, _2)]
			;

        logical =
            equality[_val = _1] >> *(logical_op > equality)[_val = ast::make_binary(_1, _val, _2)]
            ;

        equality =
            relational[_val = _1] >> -(equality_op > relational)[_val = ast::make_binary(_1, _val, _2)]
            ;

        relational =
            additive[_val = _1] >> -(relational_op > additive)[_val = ast::make_binary(_1, _val, _2)]
            ;

        additive =
            multiplicative[_val = _1] >> *(additive_op > multiplicative)[_val = ast::make_binary(_1, _val, _2)]
            ;

        multiplicative =
            factor[_val = _1] >> *(multiplicative_op > factor)[_val = ast::make_binary(_1, _val, _2)]
            ;

        factor =
            primary[_val = _1] >> *( power > factor )[_val = ast::make_binary(_1, _val, _2)]
            ;

        unary =
            (ufunc > '(' > logical > ')')[_val = ast::make_unary(_1, _2)]
            ;

        binary =
            (bfunc > '(' > logical > ',' > logical > ')')[_val = ast::make_binary(_1, _2, _3)]
            ;

        variable =
            raw[lexeme[alpha >> *(alnum | '_')]]
            ;

        primary =
              double_
            | ('(' > logical > ')')
            | (unary_op > primary)[_val = ast::make_unary(_1, _2)]
            | binary
            | unary
            | constant
            | variable
            ;

        // TODO:
        // Problem that we allow logical expressions even where
        // only double expr can appear
        // we should probably have two different primary expressions
        // and possibly even more for arrays
        // However both doouble and bool values caan be assigned to a variable
        // so we have to allow 'logical' everywhere where we allow variable
        // Thus we can really treat bool values as special double values.

        // clang-format on

        program.name("program");
        assignment.name("assignment");
        logical.name("logical");
        equality.name("equality");
        relational.name("relational");
        additive.name("additive");
        multiplicative.name("multiplicative");
        factor.name("factor");
        variable.name("variable");
        primary.name("primary");
        unary.name("unary");
        binary.name("binary");

        // typedef boost::phoenix::function<error_handler<Iterator> >
        // error_handler_function; qi::on_error<qi::fail>(expression,
        //        error_handler_function(error_handler<Iterator>())(
        //            "Error! Expecting ", qi::_4, qi::_3));
        on_error<fail>(
            program,
            boost::phoenix::bind(boost::phoenix::ref(err_handler), _3, _2, _4));

    }
};

} // namespace parser

typedef parser::grammar<std::string::const_iterator> grammar;



} // namespace bparser

#endif // INCLUDE_GRAMMAR_HH_
