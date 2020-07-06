/*
 * grammar.impl.hh
 *
 *  Created on: Feb 3, 2020
 *      Author: jb
 */

#ifndef INCLUDE_GRAMMAR_IMPL_HH_
#define INCLUDE_GRAMMAR_IMPL_HH_

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

#include "processor.hh"
#include "array.hh"
#include "ast.hh"


namespace bparser {


namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


namespace parser {




template<class T>
ast::unary_fn::function_type unary_array() {
	return static_cast<ast::unary_fn::function_type>(&(Array::unary_op<T>));
}

template<class T>
ast::binary_fn::function_type binary_array() {
	return static_cast<ast::binary_fn::function_type>(&(Array::binary_op<T>));
}

ast::binary_fn::function_type append_to() {
	return static_cast<ast::binary_fn::function_type>(&(Array::append_to));
}

ast::binary_fn::function_type none_const() {
	return static_cast<ast::binary_fn::function_type>(&(Array::append_to));
}


struct expectation_handler {
    template <typename>
    struct result {
        typedef void type;
    };

    template <typename Iterator>
    void operator()(Iterator first, Iterator last,
                    boost::spirit::info const &info) const {
    	std::ostringstream out;
    	out << info;
        Throw() << "Expected " << out.str() << " at \"" << std::string(first, last)
        << "\""; // NOLINT
    }
};

#define UN_FN(repr, fn) (repr, {repr, fn})
#define BN_FN(repr, fn) (repr, {repr, fn})

template <typename Iterator>
struct grammar : qi::grammar<Iterator, ast::operand(), ascii::space_type> {
    expectation_handler err_handler;


    qi::rule<Iterator, ast::operand(), ascii::space_type>
		program,
		definition,
		assignment,
		primary,
		enclosure,
		expression,
		conditional_expr,
		or_test,
		and_test,
		not_test,
		not_test_op,
		comparison,
		additive_expr,
		multiplicative_expr,
		signed_optional,
		signed_expr,
		power,
        array_constr,
		array_constr_list,
		literal_double,
		literal_int,
		const_lit,
		const_idx,
		subscriptable,
		subscription,
		slicing,
		slice_item,
		proper_slice,
		index_array,
		index_array_list,
		call,
		unary_call,
		binary_call
    	;
    qi::rule<Iterator, std::string()>
    	identifier;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type, ast::unary_fn>
        ufunc, unary_op, not_op;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type, ast::binary_fn>
        bfunc,
		additive_op,
		multiplicative_op,
		and_op, or_op,
		relational_op,
		equality_op,
		power_op,
		semicol_op,
		array_constr_comma_op;

    //ast::unary_op none_const, true_const, false_const;


    grammar()
    : grammar::base_type(program, "bparser_program")
    {


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
			BN_FN("//", binary_array<_mod_>())	// floor division
			BN_FN("@", binary_array<_mod_>())	// python
            BN_FN("%", binary_array<_mod_>())
            ;

        and_op.add BN_FN("and", binary_array<_and_>());
        or_op.add BN_FN("or", binary_array<_or_>());
        not_op.add UN_FN("not", unary_array<_neg_>());

        ast::binary_fn append_to_fn = {",", append_to()};
        ast::binary_fn subscribe_fn = {"[]", &Array::subscribe};
        array_constr_comma_op.add(",", append_to_fn);

//        array_slice_comma_op.add
//			BN_FN(",", &append_slice_list);

        relational_op.add
            BN_FN("<" , binary_array<_lt_>())
            BN_FN("<=", binary_array<_le_>())
            BN_FN(">" , &gt_op)
            BN_FN(">=", &ge_op)
            BN_FN("==", binary_array<_eq_>())
            BN_FN("!=", binary_array<_ne_>())
            ;

        power_op.add
            BN_FN("**", binary_array<_pow_>())
            ;



        ast::binary_fn semicol_fn = {";", &ast::semicol_fn};
        ast::ternary_fn if_else_fn = {"ifelse", &Array::if_else};
        ast::ternary_fn slice_fn = {"slice", &Array::slice};
        //auto head_const = ast::make_const("Head", &Array::empty_array);
        auto none_const = ast::make_const("None", &Array::none_array);
        auto true_const = ast::make_const("True", &Array::true_array);
        auto false_const = ast::make_const("False", &Array::false_array);



//        semicol_op.add
//            (";", &ast::semicol_fn)  // TODO: possibly some special function
//			;

        /**
         * Spirit:Qi syntax notes:
         * lhs >> rhs : match if rhs follows lhs
         * lhs > rhs  : same, but throw if there is lhs and no rhs
         * -(expr)     : optional rule
         * *(expr)     : match zero or more repetitions
         * |           : alternative (tried sequantionaly)
         *
         * Every rule have associated type it produce in  AST given as
         * second template parameter to qi::rule.
         *
         * !! If the attribute operation in bracket is specified for the part of the rule it must by specified for all its terms.
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
#define RULE(r_name) r_name.name(#r_name); r_name

        /**
         * BParser grammar tries to follow Python 3.8 grammar.
         * In order to perform all array operations statically before compilation
         * we restrict indexing and slicing to int literals instead of general expressions.
         */

        RULE(program) = (definition >> ';' >> expression)[qi::_val = ast::make_binary(semicol_fn, qi::_1,  qi::_2)]
						| expression[qi::_val = qi::_1];

        RULE(definition) = assignment[qi::_val = qi::_1] >> *(';' >> assignment)
        		           [qi::_val = ast::make_binary(semicol_fn, qi::_val, qi::_1)];

        RULE(assignment) = (identifier >> '=' > expression)[qi::_val = ast::make_assign(qi::_1, qi::_2)];








        RULE(expression) = conditional_expr.alias();
        RULE(conditional_expr) = or_test[qi::_val = qi::_1] >>
        		-(qi::lit("if") > or_test > qi::lit("else") > expression)
        		[qi::_val = ast::make_ternary(if_else_fn, qi::_val, qi::_1, qi::_2)];
        // conditional expr can only be used in the 'else' branch !!
        RULE(or_test)  = and_test[qi::_val = qi::_1] >>
        				*(or_op > and_test)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];
        RULE(and_test) = not_test[qi::_val = qi::_1] >>
        				*(and_op > not_test)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];
		RULE(not_test) = comparison | not_test_op;
		RULE(not_test_op) = (not_op > not_test)[qi::_val = ast::make_unary(qi::_1, qi::_2)];

        RULE(comparison) = additive_expr[qi::_val = qi::_1] >>
        		  	    -(relational_op > additive_expr)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];
        RULE(additive_expr) = multiplicative_expr[qi::_val = qi::_1] >>
        				*(additive_op > multiplicative_expr)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(multiplicative_expr) = signed_optional[qi::_val = qi::_1] >>
        				*(multiplicative_op > signed_optional)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(signed_optional) = signed_expr | power;
        RULE(signed_expr) = (unary_op > signed_optional)[qi::_val = ast::make_unary(qi::_1, qi::_2)];
        RULE(power) = primary[qi::_val = qi::_1] >>
        				-(power_op > signed_optional)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(primary) = literal_double | const_lit | subscription;
        RULE(subscriptable) = array_constr | enclosure | call | identifier;

        RULE(call) = binary_call | unary_call;
        RULE(unary_call) = (qi::no_skip[ufunc > '('] > expression > ')')
        		           [qi::_val = ast::make_unary(qi::_1, qi::_2)];
        RULE(binary_call) = (qi::no_skip[bfunc > '('] > expression > ',' > expression > ')')
        		            [qi::_val = ast::make_binary(qi::_1, qi::_2, qi::_3)];

        RULE(subscription) = subscriptable[qi::_val = qi::_1] >> -('[' >> slicing > ']')
                             [qi::_val = ast::make_binary(subscribe_fn, qi::_val, qi::_1)];
        RULE(slicing) = slice_item[qi::_val = ast::make_binary(append_to_fn, none_const, qi::_1)]
					    >> *(array_constr_comma_op > slice_item)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];;
        RULE(slice_item) =  index_array | proper_slice | const_idx;

        RULE(const_lit) =
        		qi::lit("None")[qi::_val = none_const] |
				qi::lit("True")[qi::_val = true_const] |
				qi::lit("False")[qi::_val = false_const];

        RULE(enclosure) = '(' > expression > ')';
        RULE(array_constr) = '[' > array_constr_list > ']';
        RULE(array_constr_list) = expression[qi::_val = ast::make_binary(append_to_fn, none_const, qi::_1)]
						 >> *(array_constr_comma_op > expression)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(literal_double) = qi::double_; // includes ints
        RULE(identifier) = qi::raw[qi::lexeme[(qi::alpha | '_') >> *(qi::alnum | '_')]];


        RULE(proper_slice) = (-(literal_int) >> qi::lit(':')
        						>> -(literal_int)
								>> -( qi::lit(':') >> literal_int))
								[qi::_val = ast::make_ternary(slice_fn,
										ast::treat_optional(qi::_1),
										ast::treat_optional(qi::_2),
										ast::treat_optional(qi::_3))];
        RULE(index_array) = '[' > index_array_list > ']';
        RULE(index_array_list) = const_idx[qi::_val = ast::make_binary(append_to_fn, none_const, qi::_1)]
						 >> *(array_constr_comma_op > const_idx)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];
        RULE(const_idx) = literal_int |  const_lit;
        RULE(literal_int) = qi::int_;



        // TODO:
        // Problem that we allow logical expressions even where
        // only double expr can appear
        // we should probably have two different primary expressions
        // and possibly even more for arrays
        // However both doouble and bool values caan be assigned to a variable
        // so we have to allow 'logical' everywhere where we allow variable
        // Thus we can really treat bool values as special double values.

        // clang-format on


        qi::on_error<qi::fail>(
            program,
            boost::phoenix::bind(boost::phoenix::ref(err_handler), qi::_3, qi::_2, qi::_4));

//        qi::debug(program);
//        qi::debug(expression);
//        qi::debug(conditional_expr);
//        qi::debug(or_test);
//        qi::debug(and_test);
//        qi::debug(not_test);
//        qi::debug(comparison);
//        qi::debug(signed_optional);
//        qi::debug(power);
//
//        qi::debug(primary);
//        qi::debug(subscription);
//        qi::debug(call);
//        qi::debug(unary_call);
//        qi::debug(binary_call);
//        qi::debug(enclosure);
//        qi::debug(primary);
//        qi::debug(identifier);
//        qi::debug(literal_double);
//        qi::debug(const_lit);
//        qi::debug(slicing);
//        qi::debug(slice_item);
//        qi::debug(const_idx);



    }
};

} // namespace parser

typedef parser::grammar<std::string::const_iterator> grammar;


} // namespace bparser

#endif /* INCLUDE_GRAMMAR_IMPL_HH_ */
