/*
 * array_ast_interface.hh
 *
 * Contains functions that are not vital part of Array direct interface, but are necessary
 * for transformation of the AST.
 *
 *  Created on: Jul 6, 2020
 *      Author: jb
 */

#ifndef INCLUDE_ARRAY_AST_INTERFACE_HH_
#define INCLUDE_ARRAY_AST_INTERFACE_HH_

#include <vector>
#include <typeinfo>
#include "array.hh"


namespace bparser {


typedef std::vector<Array> ListArray;
typedef std::vector<int> ListInt;
typedef boost::variant<int, ListInt, Slice> Range;
typedef std::vector<Range> ListRange;
typedef std::pair<Array, ListRange> SubscribePair;

/// pair( 'a<b') , b); support for  comparison chaining for the 'b'
typedef std::pair<Array, Array> ComparisonPair;


typedef boost::variant<
		int,
		double,
		Array,
		Range,
		ListRange,
		ComparisonPair
		> ParserResult;

typedef std::vector<ParserResult> ResultList;


typedef Array (*ArrayFnUnary)(const Array &);
typedef Array (*ArrayFnBinary)(const Array &, const Array &);
typedef Array (*ArrayFnTernary)(const Array &, const Array &, const Array &);
typedef Array (*ArrayFnVariadic)(const ListArray &);
typedef Array (*ArrayFnSubscribe)(const Array &, const ListRange &);
typedef Array (*ArrayFnInt)(int);
typedef Array (*ArrayFnShape)(const Shape &);
typedef Array (*ArrayFnShapeDbl)(const Shape &, double);

typedef Range (*RangeFnListInt)(const ListInt &);
typedef ListRange (*ListRangeFnListRange)(const ListRange &);
typedef Range (*RangeFnInt)(int);

struct ChainedCompareFn {
	ChainedCompareFn(ArrayFnBinary cmp_fn)
	: cmp_fn_(cmp_fn)
	{}

	ComparisonPair operator() (const ComparisonPair &pair, const Array &other) {
		const Array &previous_cmp = pair.first;
		Array current_cmp(cmp_fn_(pair.second, other));
		Array previous_and_current(Array::binary_op<details::_and_>(previous_cmp, current_cmp));
		return ComparisonPair(previous_and_current, other);
	}

	ArrayFnBinary cmp_fn_;
};

typedef Array (*CloseChain)(const ComparisonPair &);

typedef boost::variant <
	ArrayFnUnary,
	ArrayFnBinary,
	ChainedCompareFn,
	CloseChain,
	ArrayFnTernary,
	ArrayFnVariadic,
	ArrayFnSubscribe,  // slicing
	ArrayFnShape,
	ArrayFnInt,
	ArrayFnShapeDbl,
	RangeFnListInt,
	ListRangeFnListRange,  //IndexArray and Slice
	RangeFnInt
	> ArrayFn;


template <class Target, class Source>
Target get_type(Source result) {
	if (Target* l = boost::get<Target>(&result)) return *l;

	Throw() << "Internal error.\n"
		<< "Wrong type: " << typeid(Target).name() << "variant subtype:" << result.which() << "\n"
		<< "Expected: " << typeid(Target).name() << "\n";
}


template <class Target, class Source>
bool check_type(Source result) {
	if (Target* l = boost::get<Target>(&result)) return true;
	else return false;
}



/**
 * Usage: call_visitor(a_list)(fn)
 *
 * Apply function 'fn' to 'a_list' of arguments
 */
struct array_conversion_visitor {
	typedef Array result_type;
    explicit array_conversion_visitor() {}

	result_type operator()(int val) const {
    	return Array::constant({double(val)});

    }

	result_type operator()(double val) const {
    	return Array::constant({double(val)});

    }

	result_type operator()(Range val) const {
		ListInt l = get_type<ListInt>(val);
		std::vector<double> conv(l.begin(), l.end());
		Shape res_shape({uint(l.size())});
    	return Array::constant(conv, res_shape);

    }

	result_type operator()(Array val) const {
    	return val;

    }

	static result_type error(ParserResult val) {
    	Throw() << "Internal error.\n"
    		<< "Wrong type: " << "ParserResult:" << val.which() << "\n"
			<< "Expected: Array,  int, double, Shape\n";
		return Array();
	}

    result_type operator()(ListRange val) const {
		return error(val);
    }
    result_type operator()(ComparisonPair val) const {
		return error(val);
    }
};




/// Similar to get_type<Array> but with conversion from int.
static Array get_array(ParserResult result) {
	return boost::apply_visitor(array_conversion_visitor(), result);
}





/**
 * Usage: call_visitor(a_list)(fn)
 *
 * Apply function 'fn' to 'a_list' of arguments
 */
struct call_visitor {
	typedef ParserResult result_type;
	const ResultList &alist_;
    explicit call_visitor(const ResultList &alist)
	: alist_(alist)
	{}

    static ComparisonPair make_cmp_pair(ParserResult result) {
    	if (ComparisonPair* l = boost::get<ComparisonPair>(&result))
    		return *l;
    	else {
    		Array arr = get_array(result);
    		return ComparisonPair(Array::true_array(arr), arr);
    	}
    }

    result_type operator()(ArrayFnUnary fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	//std::cout << "unary fn" << "\n";
    	return fn(get_array(alist_[0]));
    }

    result_type operator()(ArrayFnBinary fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "binary fn" << "\n";
    	return fn(get_array(alist_[0]), get_array(alist_[1]));
    }

    result_type operator()(ChainedCompareFn fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "binary fn" << "\n";
    	return fn(
    			make_cmp_pair(alist_[0]),
				get_array(alist_[1]));
    }

    result_type operator()(CloseChain fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	return fn(get_type<ComparisonPair>(alist_[0]));
    }

    result_type operator()(ArrayFnTernary fn) const {
    	BP_ASSERT(alist_.size() == 3);
    	return fn(
    			get_array(alist_[0]),
				get_array(alist_[1]),
    			get_array(alist_[2]));
    }

    result_type operator()(ArrayFnVariadic fn) const {
    	std::vector<Array> alist;
    	for(ParserResult item : alist_)
    		alist.push_back(get_array(item));
    	return fn(alist);
    }

    result_type operator()(ArrayFnSubscribe fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "subscribe fn" << "\n";
    	return fn(
    			get_type<Array>(alist_[0]),
				get_type<ListRange>(alist_[1]));
    }

    result_type operator()(ArrayFnShape fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	//std::cout << "int list fn" << "\n";
    	Range range =get_type<Range>(alist_[0]);
		ListInt l = get_type<ListInt>(range);
    	Shape shape(l.begin(), l.end());
    	return fn(shape);
    }

    result_type operator()(ArrayFnShapeDbl fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	Range range = get_type<Range>(alist_[0]);
		ListInt l = get_type<ListInt>(range);
    	Shape shape(l.begin(), l.end());
    	return fn(shape, get_type<double>(alist_[1]));
    }

    result_type operator()(ArrayFnInt fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	return fn(get_type<int>(alist_[0]));
    }

    result_type operator()(RangeFnListInt fn) const {
    	std::vector<int> alist;
    	//std::cout << "int list fn" << "\n";
    	for(ParserResult item : alist_)
    		alist.push_back(get_type<int>(item));
    	return fn(alist);
    }

    result_type operator()(ListRangeFnListRange fn) const {
    	// convert ListInt to ListRange
    	ListRange r_list;
    	//std::cout << "range list fn" << "\n";
    	for(auto item : alist_) {
    		r_list.push_back(get_type<Range>(item));
    	}

    	return fn(r_list);
    }

    result_type operator()(RangeFnInt fn) const {
    	//std::cout << "index fn" << "\n";
    	return fn(get_type<int>(alist_[0]));
    }

};


struct subscribe_visitor {
	typedef MultiIdxRange result_type;
	const MultiIdxRange &range_;
	uint axis_;
	subscribe_visitor(const MultiIdxRange &r, uint axis)
	: range_(r), axis_(axis)
	{}

	result_type operator()(int idx) const {
		MultiIdxRange res(range_);
		res.sub_index(axis_, idx);
		return res;
	}

	result_type operator()(Slice slice_range) const {
		MultiIdxRange res(range_);
		res.sub_slice(axis_, slice_range);
		return res;
	}

	result_type operator()(ListInt index_range) const {
		MultiIdxRange res(range_);
		res.sub_range(axis_, index_range);
		return res;
	}

};

struct NamedArrayFn {
	std::string repr;
	ArrayFn fn;
};

inline Range empty_array(const ListInt& UNUSED(x)) {
	Throw() << "Empty Array not allowed.";
	return ListInt();
}


inline Range create_index_array(const ListInt& l) {
	return l;
}

inline Range create_slice(const ListInt& l) {
	Slice s;
	std::copy_n(std::make_move_iterator(l.begin()), 3, s.begin());
	return s;
}

inline Range create_index(int idx) {
	return idx;
}

inline Array subscribe(const Array &a, const ListRange &slice_list) {
	//std::cout << "subscribe start" << slice_list.size() << "\n";
	auto mir = MultiIdxRange(a.shape()).full();
	for(uint i=0; i<slice_list.size(); i++) {
		Range axis_range = slice_list[i];
		mir = boost::apply_visitor(subscribe_visitor(mir, i), axis_range);
	}
	//std::cout << "sub shape:  " << print_vector(mir.sub_shape()) << "\n";
	//std::cout << "full shape: " << print_vector(mir.full_shape_) << "\n";
	//std::cout << "ranges: " << print_vector(mir.ranges_[0]) << "\n";

	return Array(a, mir);  // TODO: make valid implementation
}

inline ListRange range_list(const ListRange &slice_list) {
	return slice_list;
}

inline Array close_chain(const ComparisonPair & x) {
	return x.first;
}


} // namespace bparser
#endif /* INCLUDE_ARRAY_AST_INTERFACE_HH_ */
