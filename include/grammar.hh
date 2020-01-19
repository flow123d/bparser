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

#include "array.hh"
#include "ast.hh"


namespace bparser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


namespace parser {




template<class T>
ast::unary_fn::type unary_array() {
	return static_cast<ast::unary_fn::type>(&(Array::unary_op<T>));
}

template<class T>
ast::binary_fn::type binary_array() {
	return static_cast<ast::binary_fn::type>(&(Array::binary_op<T>));
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

        Throw(msg.str()); // NOLINT
    }
};

#define UN_FN(repr, fn) (repr, {repr, fn})
#define BN_FN(repr, fn) (repr, {repr, fn})

template <typename Iterator>
struct grammar : qi::grammar<Iterator, ast::operand(), ascii::space_type> {
    expectation_handler err_handler;


    qi::rule<Iterator, ast::operand(), ascii::space_type>
    	primary, logical, equality, relational, additive, multiplicative, factor, unary, binary, program, definition, assignment, unary_op_r;
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
            UN_FN("abs"  , unary_array<_abs_>())
        	UN_FN("acos" , unary_array<_acos_>())
			UN_FN("asin" , unary_array<_asin_>())
			UN_FN("atan" , unary_array<_atan_>())
			UN_FN("ceil" , unary_array<_ceil_>())
            UN_FN("cos"  , unary_array<_cos_>())
            UN_FN("cosh" , unary_array<_cosh_>())
            UN_FN("deg"  , &deg_fn)
            UN_FN("exp"  , unary_array<_exp_>())
            UN_FN("floor", unary_array<_floor_>())
            UN_FN("isinf", unary_array<_isinf_>())
            UN_FN("isnan", unary_array<_isnan_>())
            UN_FN("log"  , unary_array<_log_>())
            UN_FN("log10", unary_array<_log10_>())
            UN_FN("rad"  , &rad_fn)
            UN_FN("sgn"  , unary_array<_sgn_>())
            UN_FN("sin"  , unary_array<_sin_>())
            UN_FN("sinh" , unary_array<_sinh_>())
            UN_FN("sqrt" , unary_array<_sqrt_>())
            UN_FN("tan"  , unary_array<_tan_>())
            UN_FN("tanh" , unary_array<_tanh_>())
            ;

        bfunc.add
            BN_FN("atan2", binary_array<_atan2_>())
            BN_FN("pow"  , binary_array<_pow_>())
            ;

        unary_op.add
            UN_FN("+", &unary_plus)
            UN_FN("-", unary_array<_minus_>())
            UN_FN("!", unary_array<_neg_>())
            ;

        additive_op.add
            BN_FN("+", binary_array<_add_>())
            BN_FN("-", binary_array<_sub_>())
            ;

        multiplicative_op.add
            BN_FN("*", binary_array<_mul_>())
            BN_FN("/", binary_array<_div_>())
            BN_FN("%", binary_array<_mod_>())
            ;

        logical_op.add
            BN_FN("&&", binary_array<_and_>())
            BN_FN("||", binary_array<_or_>())
            ;

        relational_op.add
            BN_FN("<" , binary_array<_lt_>())
            BN_FN("<=", binary_array<_le_>())
            BN_FN(">" , &gt_op)
            BN_FN(">=", &ge_op)
            BN_FN("==", binary_array<_eq_>())
            BN_FN("!=", binary_array<_ne_>())
            ;


        power.add
            BN_FN("**", binary_array<_pow_>())
            ;

        ast::binary_fn semicol_fn = {";", &ast::semicol_fn};

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
        program =
        	(definition >> ';' >> logical)[_val = ast::make_binary(semicol_fn, _1,  _2)]
			| logical[_val = _1]
        //program = -(definition > ';')[_val = _1] >> logical[_val = _1]
			;




        //program = definition.alias();
        definition =
        		assignment[_val = _1] >> *(';' >> assignment)[_val = ast::make_binary(semicol_fn, _val, _1)]
				;

        assignment =
        	(variable >> '=' > logical)[_val = ast::make_assign(_1, _2)]
			;

        //program = logical.alias();

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
            primary[_val = _1] >> -(power > factor)[_val = ast::make_binary(_1, _val, _2)]

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

        unary_op_r = (unary_op > primary)[_val = ast::make_unary(_1, _2)];

        primary =
              double_
            | ('(' > logical > ')')
            | unary_op_r
            | binary
            | unary
            | constant
            | variable
			| eps
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

        on_error<fail>(
            program,
            boost::phoenix::bind(boost::phoenix::ref(err_handler), _3, _2, _4));

    }
};

} // namespace parser

typedef parser::grammar<std::string::const_iterator> grammar;



} // namespace bparser

#endif // INCLUDE_GRAMMAR_HH_
