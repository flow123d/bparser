/*
 * parser.hh
 *
 *  Created on: Dec 26, 2019
 *      Author: jb
 */

#ifndef INCLUDE_PARSER_HH_
#define INCLUDE_PARSER_HH_

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_utree.hpp>


#include <iostream>
#include <string>

namespace bparser
{
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	namespace spirit = boost::spirit;




	/**
	 * Example:
	 * SymbolTable st;
	 * st.add("pressure", pressure_array)
	 */
	struct SymbolTable : qi::symbols<char, Array> {

	};



    ///////////////////////////////////////////////////////////////////////////////
    //  Our calculator grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct grammar : qi::grammar<Iterator, ascii::space_type, spirit::utree()>
    {
        grammar() : grammar::base_type(expression)
        {
            using qi::uint_;
            using qi::char_;

            Rule mult_operand =
                    double_
                    |   '(' >> add_expr >> ')';

            Rule mult_expr =
                    mult_first_opearand
                    |   mult_op
            		|   div_op;

            Rule mult_op =
                mult_operand >> char_('*') >> mult_expr;

            Rule div_op =
                mult_first_operand >> char_('/') >> mult_expr;

            Rule add_first_operand =
            		mult_expr
					| sign_plus
					| sign_minus;

            Rule sign_plus =
            		char_('+') >> mult_expr;

            Rule sign_minus =
            		char_('-') >> mult_expr;


       		Rule add_expr =
            		mult_op
					| add_op
            		| sub_op;

            Rule add_op =
            		add_first_operand >> char_('+') >> add_expr;

            Rule sub_op =
            		add_first_operand >> char_('-') >> add_expr;


            expression = add_expr;
//            BOOST_SPIRIT_DEBUG_NODE(expression);
//            BOOST_SPIRIT_DEBUG_NODE(term);
//            BOOST_SPIRIT_DEBUG_NODE(factor);
        }

        qi::rule<Iterator, ascii::space_type, spirit::utree()> expression;
    };
}



#endif /* INCLUDE_PARSER_HH_ */
