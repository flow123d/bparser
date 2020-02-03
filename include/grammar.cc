/*
 * grammar.cc
 *
 *  Created on: Feb 3, 2020
 *      Author: jb
 */

#include <string>
#include "ast.hh"
#include "grammar.hh"
#include "grammar.impl.hh"

namespace bparser {

bool parse_expr(std::string expr, ast::operand &ast) {
    std::string::const_iterator first = expr.begin();
    std::string::const_iterator last = expr.end();

    boost::spirit::ascii::space_type space;
    bool r = qi::phrase_parse(
        first, last, grammar(), space,
        ast);

    if (!r || first != last) {
        std::string rest(first, last);
        Throw("Parsing failed at " + rest); // NOLINT
    }

}


} // namespace bparser
