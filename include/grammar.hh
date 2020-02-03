#ifndef INCLUDE_GRAMMAR_HH_
#define INCLUDE_GRAMMAR_HH_


#include <string>
#include "ast.hh"


namespace bparser {

bool parse_expr(std::string expr, ast::operand &ast);

} // namespace bparser

#endif // INCLUDE_GRAMMAR_HH_
