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
//#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix.hpp>


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
ArrayFn unary_array() {
	return static_cast<ArrayFnUnary>(&(Array::unary_op<T>));
}

template<class T>
ArrayFn binary_array() {
	return static_cast<ArrayFnBinary>(&(Array::binary_op<T>));
}

Array error_reserved(const Array & UNUSED(x)) {
	// TODO report
	Throw() << "Reserved identifier.";
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

#define FN(repr, fn) (repr, {repr, fn})

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
		comparison_chain,
		comparison_chained,
		additive_expr,
		multiplicative_expr,
		signed_optional,
		signed_expr,
		power,
        array_constr,
		array_constr_list_opt,
		array_constr_list,
		literal_number,
		literal_int,
		const_lit,
		const_idx,
		const_index,
		subscriptable,
		subscription,
		slicing,
		slice_item,
		proper_slice,
		index_array,
		index_array_list,
		call,
		param_list_opt,
		param_list,
		param,
		unary_call,
		binary_call,
		none_lit
    	;
    qi::rule<Iterator, std::string()>
    	identifier;

    qi::symbols<typename std::iterator_traits<Iterator>::value_type, NamedArrayFn>
    	reserved,
		func,
		unary_op,
		not_op,
		additive_op,
		multiplicative_op,
		and_op, or_op,
		relational_op,
		equality_op,
		power_op,
		semicol_op;
    //ast::unary_op none_const, true_const, false_const;


    grammar()
    : grammar::base_type(program, "bparser_program")
    {
    	reserved.add
    			FN("None", &error_reserved)
    			;

        func.add
            FN("abs"  , unary_array<_abs_>())
        	FN("acos" , unary_array<_acos_>())
			FN("asin" , unary_array<_asin_>())
			FN("atan" , unary_array<_atan_>())
			FN("ceil" , unary_array<_ceil_>())
            FN("cos"  , unary_array<_cos_>())
            FN("cosh" , unary_array<_cosh_>())
            FN("rad2deg"  , &deg_fn)
            FN("exp"  , unary_array<_exp_>())
            FN("floor", unary_array<_floor_>())
            FN("isinf", unary_array<_isinf_>())	// possibly replaced by inf constant
            FN("isnan", unary_array<_isnan_>()) // possibly replaced by nan constant
            FN("log"  , unary_array<_log_>())
            FN("log10", unary_array<_log10_>())
            FN("log2", unary_array<_log2_>())
            FN("deg2rad"  , &rad_fn)
            FN("sgn"  , unary_array<_sgn_>())
            FN("sin"  , unary_array<_sin_>())
            FN("sinh" , unary_array<_sinh_>())
            FN("sqrt" , unary_array<_sqrt_>())
            FN("tan"  , unary_array<_tan_>())
            FN("tanh" , unary_array<_tanh_>())
			FN("flatten", &Array::flatten)
			FN("eye"  , &Array::eye)
			FN("zeros"  , &Array::zeros)
			FN("ones"  , &Array::ones)
			FN("full"  , &Array::full)
            FN("atan2", binary_array<_atan2_>())
            FN("power"  , binary_array<_pow_>())
			FN("minimum", binary_array<_min_>())
			FN("maximum", binary_array<_max_>())
            ;

        unary_op.add
            FN("+", &unary_plus)
            FN("-", unary_array<_minus_>())
            ;

        additive_op.add
            FN("+", binary_array<_add_>())
            FN("-", binary_array<_sub_>())
            ;

        multiplicative_op.add
            FN("*", binary_array<_mul_>())
            FN("/", binary_array<_div_>())
			FN("//", &floor_div)	// floor division
			FN("@", &Array::mat_mult)	// numpy matrix multiplication
            FN("%", binary_array<_mod_>())
            ;

        and_op.add FN("and", binary_array<_and_>());
        or_op.add FN("or", binary_array<_or_>());
        not_op.add FN("not", unary_array<_neg_>());

        NamedArrayFn subscribe_fn = {"[]", &subscribe};

        relational_op.add
            FN("<" , ChainedCompareFn(Array::binary_op<_lt_>))
            FN("<=", ChainedCompareFn(Array::binary_op<_le_>))
            FN(">" , ChainedCompareFn(&gt_op))
            FN(">=", ChainedCompareFn(&ge_op))
            FN("==", ChainedCompareFn(Array::binary_op<_eq_>))
            FN("!=", ChainedCompareFn(Array::binary_op<_ne_>))
            ;

        power_op.add
            FN("**", binary_array<_pow_>())
            ;


        //NamedArrayFn array_fn = {"array", &Array::stack_zero};
        NamedArrayFn empty_array_fn = {"empty_array", &empty_array};
        NamedArrayFn semicol_fn = {";", &ast::semicol_fn};
        NamedArrayFn if_else_fn = {"ifelse", &Array::if_else};
        NamedArrayFn slice_fn = {"slice", &create_slice};
        NamedArrayFn index_fn = {"index", &create_index};
        NamedArrayFn range_list_fn = {"list", &range_list};
        NamedArrayFn close_chain_fn = {"close_chain", &close_chain};

        auto true_const = ast::make_const("True", &Array::true_array);
        auto false_const = ast::make_const("False", &Array::false_array);

        qi::real_parser<double, qi::strict_real_policies<double>> strict_double;

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
		RULE(not_test) = not_test_op | comparison_chain;
		RULE(not_test_op) = (not_op > not_test)[qi::_val = ast::make_unary(qi::_1, qi::_2)];

 		RULE(comparison_chain) =  comparison_chained | additive_expr;
        RULE(comparison_chained) = (additive_expr[qi::_val = qi::_1]
        		     >> +(relational_op > additive_expr)[qi::_val = ast::make_chained_comparison(qi::_1, qi::_val, qi::_2)]
						)[qi::_val = ast::make_unary( close_chain_fn ,qi::_val)];

        RULE(additive_expr) = multiplicative_expr[qi::_val = qi::_1] >>
        				*(additive_op > multiplicative_expr)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];
        RULE(multiplicative_expr) = signed_optional[qi::_val = qi::_1] >>
        				*(multiplicative_op > signed_optional)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(signed_optional) = signed_expr | power;
        RULE(signed_expr) = (unary_op > signed_optional)[qi::_val = ast::make_unary(qi::_1, qi::_2)];
        RULE(power) = primary[qi::_val = qi::_1] >>
        				-(power_op > signed_optional)[qi::_val = ast::make_binary(qi::_1, qi::_val, qi::_2)];

        RULE(primary) = literal_number | const_lit | subscription;
        RULE(subscriptable) = array_constr | enclosure | call | identifier;

        RULE(call) = (qi::no_skip[func > '('] > param_list_opt > ')')
        		           [qi::_val = ast::make_call(qi::_1, qi::_2)];
        RULE(param_list_opt) = (-(param_list))
        		[qi::_val = ast::treat_optional(qi::_1, none_int)]; //TODO: test fail for function with no arguments.
        RULE(param_list) = param[qi::_val = ast::make_list(0.0, qi::_1)]
						 >> *(',' > param)[qi::_val = ast::make_list(qi::_val, qi::_1)];

        RULE(param) = expression;
        RULE(subscription) = subscriptable[qi::_val = qi::_1] >> -('[' >> slicing > ']')
                             [qi::_val = ast::make_binary(subscribe_fn, qi::_val,
                            		 ast::make_call(range_list_fn, qi::_1))];
        RULE(slicing) = slice_item[qi::_val = ast::make_list(0.0, qi::_1)]
					    >> *("," > slice_item)[qi::_val = ast::make_list(qi::_val, qi::_1)];;
        RULE(slice_item) =  index_array | proper_slice | const_index;

        RULE(const_lit) =
				qi::lit("True")[qi::_val = true_const] |
				qi::lit("False")[qi::_val = false_const];

        RULE(enclosure) = '(' > expression > ')';
        RULE(array_constr) = ('[' > array_constr_list_opt  > ']')
        		             [qi::_val = qi::_1];
        RULE(array_constr_list_opt) = (-(array_constr_list))
        		[qi::_val = ast::treat_optional(qi::_1, ast::make_unary(empty_array_fn, none_int))];
        RULE(array_constr_list) = expression[qi::_val = ast::make_list(0.0, qi::_1)]
						 >> *(',' > expression)[qi::_val = ast::make_list(qi::_val, qi::_1)];
        RULE(literal_number) = (strict_double | qi::int_); // includes ints

        RULE(identifier) =
        		&(! reserved) >>
        		qi::raw[qi::lexeme[(qi::alpha | '_') >> *(qi::alnum | '_')]];


        RULE(proper_slice) = (-(literal_int) >> qi::lit(':')
        						>> -(literal_int)
								>> -( qi::lit(':') >> literal_int))
								[qi::_val = ast::make_ternary(slice_fn,
										ast::treat_optional(qi::_1, none_int),
										ast::treat_optional(qi::_2, none_int),
										ast::treat_optional(qi::_3, none_int))];
        RULE(index_array) = ('[' > index_array_list > ']')
						[qi::_val = qi::_1];
        RULE(index_array_list) = const_idx[qi::_val = ast::make_list(none_int, qi::_1)]
						 >> *("," > const_idx)[qi::_val = ast::make_list(qi::_val, qi::_1)];
        RULE(const_index) = (const_idx)[qi::_val = ast::make_unary(
        		index_fn, qi::_1)];
        RULE(const_idx) = (literal_int |  none_lit);
        RULE(none_lit) = qi::lit("None")[qi::_val = none_int];
        //RULE(const_idx) = (literal_int |  const_lit)[qi::_val = ast::make_call(index_array_fn, qi::_1)];

        // TODO: test that we can return an index array with single index
        // TODO: single index array works, but we need distinction between a[1] and a[[1]].
        // Former case reduce the shape, while the later does not.
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
//        qi::debug(comparison_chained);
//        qi::debug(signed_optional);
//        qi::debug(power);
//        qi::debug(array_constr);
//        qi::debug(array_constr_list);
//        qi::debug(index_array);
//        qi::debug(index_array_list);
//
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
//        none_lit


    }
};

} // namespace parser

typedef parser::grammar<std::string::const_iterator> grammar;


} // namespace bparser

#endif /* INCLUDE_GRAMMAR_IMPL_HH_ */
