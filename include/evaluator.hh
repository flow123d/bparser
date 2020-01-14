



#include <map>
#include <set>
#include <string>
#include <vector>
#include <boost/variant/get.hpp>

#include "ast.hh"
#include "assert.hh"

namespace bparser {

namespace ast {

//// Optimizer
//
//template <typename T, typename U>
//struct is_same {
//    static const bool value = false;
//};
//
//template <typename T>
//struct is_same<T, T> {
//    static const bool value = true;
//};
//
//template <typename T>
//struct holds_alternative_impl {
//    typedef bool result_type;
//
//    template <typename U>
//    bool operator()(U const & /*unused*/) const {
//        return is_same<U, T>::value;
//    }
//};
//
//
///**
// * Returns true if operand is type T.
// * TODO: try to use a type trait to check that
// */
//template <typename T>
//bool holds_alternative(operand const &v) {
//    return boost::apply_visitor(holds_alternative_impl<T>(), v);
//}
//
//
//
//
//
//struct ConstantFolder {
//    typedef operand result_type;
//
//
//
//    result_type operator()(nil) const { return 0; }
//
//    result_type operator()(double n) const {
//        return n;
//    }
//
//
//    result_type operator()(std::string const &c) const {
//        return c;
//    }
//
//    result_type operator()(operation const &x, operand const &lhs) const {
//        operand rhs = boost::apply_visitor(*this, x.rhs);
//
//        if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
//            return x.op(boost::get<double>(lhs), boost::get<double>(rhs));
//        }
//        return binary_op(x.op, lhs, rhs);
//    }
//
//    result_type operator()(unary_op const &x) const {
//        operand rhs = boost::apply_visitor(*this, x.rhs);
//
//        /// If the operand is known, we can directly evaluate the function.
//        if (holds_alternative<double>(rhs)) {
//            return x.op(boost::get<double>(rhs));
//        }
//        return unary_op(x.op, rhs);
//    }
//
//    result_type operator()(binary_op const &x) const {
//        operand lhs = boost::apply_visitor(*this, x.lhs);
//        operand rhs = boost::apply_visitor(*this, x.rhs);
//
//        /// If both operands are known, we can directly evaluate the function,
//        /// else we just update the children with the new expressions.
//        if (holds_alternative<double>(lhs) && holds_alternative<double>(rhs)) {
//            return x.op(boost::get<double>(lhs), boost::get<double>(rhs));
//        }
//        return binary_op(x.op, lhs, rhs);
//    }
//
//    result_type operator()(expression const &x) const {
//        operand state = boost::apply_visitor(*this, x.lhs);
//        for (std::list<operation>::const_iterator it = x.rhs.begin();
//             it != x.rhs.end(); ++it) {
//            state = (*this)(*it, state);
//        }
//        return state;
//    }
//
//};
//
//
//
//

//
//
//struct make_array {
//    typedef expr::Array result_type;
//
//
//    explicit make_array()
//    {}
//
//
//    result_type operator()(nil) const {
//        BOOST_ASSERT(0);
//        return expr::Array();
//    }
//
//    result_type operator()(double x) const
//    {
//    	return expr::Array::constant({x});
//    }
//
//    result_type operator()(std::string const &x) const  {
//        auto it = expr_defs.find(x);
//        if (it == expr_defs.end()) {
//        	return it->second;
//        } else {
//        	std::ostringstream s;
//        	s << "Undefined var: " << x << "\n";
//        	// We do not call visitor for the assign_op so this must be error.
//        	Throw(s.str());
//        }
//    }
//
//    result_type operator()(unary_op const &x) const {
//        return x.op(boost::apply_visitor(*this, x.rhs));
//    }
//
//    result_type operator()(binary_op const &x) const {
//    	result_type lhs = boost::apply_visitor(*this, x.lhs);
//    	result_type rhs = boost::apply_visitor(*this, x.rhs);
//        return x.op(lhs, rhs);
//    }
//
//    result_type operator()(assign_op x) const  {
//    	//std::string& var_name = boos);
//    	result_type rhs = boost::apply_visitor(*this, x.rhs);
//    	expr_defs[x.lhs] = rhs;
//        return rhs;
//    }
//
//    result_type operator()(operation const &x, result_type lhs) const {
//        return x.op(lhs, boost::apply_visitor(*this, x.rhs));
//    }
//
//    result_type operator()(expression const &x) const {
//    	result_type state = boost::apply_visitor(*this, x.lhs);
//        for (std::list<operation>::const_iterator it = x.rhs.begin();
//             it != x.rhs.end(); ++it) {
//        	state = (*this)(*it, state);
//        }
//        return state;
//    }
//
//
//private:
//    mutable std::map<std::string, expr::Array> expr_defs;
//};
//
//
//
//
struct get_variables {
    typedef std::vector<std::string> result_type;

    static result_type merge(result_type a, result_type b) {
    	result_type res(a.begin(), a.end());
    	res.insert(res.end(), b.begin(), b.end());
    	return res;
    }

    explicit get_variables()
    {}


    result_type operator()(nil) const {
        BOOST_ASSERT(0);
        return {};
    }

    result_type operator()(double x) const
    { return {}; }

    result_type operator()(std::string const &x) const  {
        auto it = expr_defs.find(x);
        if (it == expr_defs.end()) {
        	return {x};
        }
        return {};
    }


    result_type operator()(unary_op const &x) const {
        return boost::apply_visitor(*this, x.rhs);
    }

    result_type operator()(binary_op const &x) const {
    	result_type lhs = boost::apply_visitor(*this, x.lhs);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
        return merge(lhs, rhs);
    }

    result_type operator()(assign_op const &x) const  {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	expr_defs.insert(x.lhs);
        return rhs;
    }



private:
    mutable std::set<std::string> expr_defs;
};


struct remove_nil {
    typedef operand result_type;


    explicit remove_nil()
    {}


    result_type operator()(nil) const {
        return nil();
    }

    result_type operator()(double x) const
    { return x; }

    result_type operator()(std::string const &x) const  {
        return x;
    }


    result_type operator()(unary_op const &x) const {
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	if (rhs.type() != typeid(nil)) return x;
    	return rhs;
    }

    result_type operator()(binary_op const &x) const {
    	result_type lhs = boost::apply_visitor(*this, x.lhs);
    	if (lhs.type() != typeid(nil)) return x;
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	if (rhs.type() != typeid(nil)) return x;
    	return rhs;
    }

    result_type operator()(assign_op const &x) const  {
    	ASSERT(x.lhs.size() > 0);
    	result_type rhs = boost::apply_visitor(*this, x.rhs);
    	ASSERT(rhs.type() != typeid(nil));
        return x;
    }



};


} // namespace ast

} // namespace matheval


