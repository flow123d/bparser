/*
 * grammar.hh
 *
 *  Created on: Dec 26, 2019
 *      Author: jb
 */

#ifndef INCLUDE_GRAMMAR_OLD_HH_
#define INCLUDE_GRAMMAR_OLD_HH_


#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_utree.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <iostream>
#include <string>

#define BOOST_SPIRIT_DEBUG

namespace bparser
{
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	namespace spirit = boost::spirit;
	namespace phx = boost::phoenix;

/**
 * Spirit summary
 *
 */


//	/**
//	 * Example:
//	 * SymbolTable st;
//	 * st.add("pressure", pressure_array)
//	 */
//	struct SymbolTable : qi::symbols<char, Array> {
//
//	};

	struct Functions {

	};

    ///////////////////////////////////////////////////////////////////////////////
    //  Our calculator grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct grammar : qi::grammar<Iterator, ascii::space_type>
    {
    	typedef qi::rule<Iterator, ascii::space_type> Rule;
        grammar() : grammar::base_type(expr, "expression_grammar")
        {
        	init_function_symbol();
        	init_variable_symbol();

//            mult_operand.name("mult_div_operand");
//            mult_operand =
//                    qi::double_
//                    ;//|   '(' >> add_expr >> ')';
//
////            mult_op =
////                mult_operand >> char_('*') >> mult_expr;
//            mult_op =
//                mult_expr >> char_('*') >> mult_operand;
//
////            div_op =
////                mult_operand >> char_('/') >> mult_expr;
//
//            mult_expr =
//                    mult_op
////            		|   div_op
//            		| mult_operand;
//
//
//            sign_plus =
//            		char_('+') >> mult_expr;
//
//            sign_minus =
//            		char_('-') >> mult_expr;
//
//            add_operand =
//            		mult_expr
//					| sign_plus
//					| sign_minus;
//
//            add_op =
//            		add_operand >> char_('+') >> add_expr;
//
//            sub_op =
//            		add_operand >> char_('-') >> add_expr;
//
//       		add_expr =
//            		mult_op
//					| add_op
//            		| sub_op;
//
//
//
//            //expression = add_expr;
//       		expression = mult_op;
//            BOOST_SPIRIT_DEBUG_NODE(expression);
////            BOOST_SPIRIT_DEBUG_NODE(term);
////            BOOST_SPIRIT_DEBUG_NODE(factor);
//
//            qi::on_error<qi::fail>
//             (
//                 expression
//               , std::cout
//                     << phx::val("Error! Expecting ")
//                     << qi::_4                               // what failed?
//                     << phx::val(" here: \"")
//                     << phx::construct<std::string>(qi::_3, qi::_2)   // iterators to error-pos, end
//                     << phx::val("\"")
//                     << std::endl
//             );

        	//assign = variable_symbol >> qi::char_('=') >> expr;

        	expr.name("expression");
            expr =
            		mult_op >> add_tail;

            add_tail.name("add_tail");
            add_tail =
                   (qi::char_('+') >> mult_op >> add_tail)
                  |(qi::char_('-') >> mult_op >> add_tail)
				  | qi::eps;


            mult_op.name("mult_op");
            mult_op = pow_op >> mult_tail;

            mult_tail.name("mult_tail");
            mult_tail =

					  (qi::char_('*') >> pow_op >> mult_tail)
					| (qi::char_('/') >> pow_op >> mult_tail)
					| (qi::char_('%') >> pow_op >> mult_tail)
					| (qi::char_('@') >> pow_op >> mult_tail)
					| qi::eps
					;

            pow_op.name("pow_op");
            pow_op = unary_op >> pow_tail;

            pow_tail.name("pow_tail");
            pow_tail =
            		 qi::char_('*') >> qi::char_('*') >> pow_op
					 | qi::eps
					 ;

            unary_op.name("unary_op");
            unary_op =
            		 (qi::char_('+') >> no_op)
            		| (qi::char_('-') >> no_op)
					| no_op
					;

            no_op.name("no_op");
            no_op =
            	qi::double_[action_double]
//				| fn_call
//				| variable_symbol
                |   '(' >> expr >> ')'
                ;

//            fn_call = function_symbol >> qi::char_('(') >> arg_list >> qi::char_(')');
//
//            arg_list = expr >> arg_tail;
//
//            arg_tail =
//					(qi::char_(',') >> arg_tail)
//            		| qi::eps;


            debug(expr);
            debug(add_tail);
            debug(mult_op);
            debug(mult_tail);
            debug(pow_op);
            debug(pow_tail);
            debug(unary_op);
            debug(no_op);

            //BOOST_SPIRIT_DEBUG_NODE(add_tail);
            //BOOST_SPIRIT_DEBUG_NODE(factor);
        }

        void init_function_symbol() {
        	function_symbol.add
				("abs", 0)		// single argument
				("sin", 3)
				("cos", 4)
				("sin", 5)
				("min", 1)		// two (or more arguments)
				("max", 2);

        }

        void init_variable_symbol() {
        	variable_symbol.add
				("pi", 0)		// single argument
				("e", 1);

        }

        void action_double(double x) {
        	std::cout << "Array::constant(" << x << ")";
        }

        qi::symbols<char, int> variable_symbol;
        qi::symbols<char, int> function_symbol;
        Rule expr, add_tail, mult_op, mult_tail, pow_op, pow_tail, unary_op, no_op, fn_call, arg_list, arg_tail;


    };


}



#endif /* INCLUDE_GRAMMAR_OLD_HH_ */
