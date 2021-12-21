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


typedef std::vector<Array> ArrayList;
typedef std::vector<int> IndexList;
typedef boost::variant<int, IndexList, Slice> Range;
typedef std::vector<Range> RangeList;
typedef std::pair<Array, RangeList> SubscribePair;

/// pair( 'a<b') , b); support for  comparison chaining for the 'b'
typedef std::pair<Array, Array> ComparisonPair;

//typedef boost::variant<
//		ArrayList,
//		IndexList,
//		Slice,
//		RangeList,
//		SubscribePair,
//		> ResultList;

typedef boost::variant<
		Array,
		int,
		Range,
		RangeList,
		ComparisonPair
		> ParserResult;

typedef std::vector<ParserResult> ResultList;

//struct append_visitor : boost::static_visitor<ResultList> {
//	ResultList head_;
//
//	explicit append_visitor(ResultList head)
//	: head_(head)
//	{}
//
//
//
//
//
//	Slice operator()(const Slice& alist) const {
//		BP_ASSERT(alist.size() == 1);
//		ResultList res(head_);
//		if (Slice *r = boost::get<Slice>(&res)) {
//			return *r;
//		} else {
//			BP_ASSERT(false);
//		}
//	}
//
//	Slice operator()(const Slice& alist) const {
//		BP_ASSERT(alist.size() == 1);
//		ResultList res(head_);
//		if (Slice *r = boost::get<Slice>(&res)) {
//			return *r;
//		} else {
//			BP_ASSERT(false);
//		}
//	}
//
//
//	template<class T>
//	ResultList operator()(const T& alist) const {
//		BP_ASSERT(alist.size() == 1);
//		ResultList res(head_);
//		if (T *r = boost::get<T>(&res)) {
//			r->push_back(alist[0]);
//			return res;
//		} else {
//			BP_ASSERT(false);
//		}
//	}
//};

//struct size_visitor : boost::static_visitor<uint> {
//	template<class T>
//	uint operator()(const T& alist) const {
//		return alist.size();
//	}
//
//};
//
//inline uint result_size(ResultList r) {
//	return boost::apply_visitor(size_visitor(), r);
//}





typedef Array (*ArrayFnUnary)(const Array &);
typedef Array (*ArrayFnBinary)(const Array &, const Array &);
typedef Array (*ArrayFnTernary)(const Array &, const Array &, const Array &);
typedef Array (*ArrayFnVariadic)(const ArrayList &);
typedef Range (*IntListFn)(const IndexList &);
typedef Array (*SubscribeFn)(const Array &, const RangeList &);
typedef RangeList (*RangeListFn)(const RangeList &);
typedef Range (*IndexFn)(int);

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
	IntListFn,  //IndexArray and Slice
	SubscribeFn,  // slicing
	RangeListFn,
	IndexFn
	> ArrayFn;

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

    template <class T>
    static T get_type(ParserResult result) {
    	if (T* l = boost::get<T>(&result)) return *l;

    	Throw() << "Internal error.\n"
    		<< "Wrong type: " << "ParserResult:" << result.which() << "\n"
			<< "Expected: " << typeid(T).name() << "\n";
    }

    static ComparisonPair make_cmp_pair(ParserResult result) {
    	if (ComparisonPair* l = boost::get<ComparisonPair>(&result))
    		return *l;
    	else
    		return ComparisonPair(Array::true_array(get_type<Array>(result)), get_type<Array>(result));
    }

    result_type operator()(ArrayFnUnary fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	//std::cout << "unary fn" << "\n";
    	return fn(get_type<Array>(alist_[0]));
    }

    result_type operator()(ArrayFnBinary fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "binary fn" << "\n";
    	return fn(
    			get_type<Array>(alist_[0]),
				get_type<Array>(alist_[1]));
    }

    result_type operator()(ChainedCompareFn fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "binary fn" << "\n";
    	return fn(
    			make_cmp_pair(alist_[0]),
				get_type<Array>(alist_[1]));
    }

    result_type operator()(CloseChain fn) const {
    	BP_ASSERT(alist_.size() == 1);
    	return fn(get_type<ComparisonPair>(alist_[0]));
    }

    result_type operator()(ArrayFnTernary fn) const {
    	BP_ASSERT(alist_.size() == 3);
    	return fn(
    			get_type<Array>(alist_[0]),
				get_type<Array>(alist_[1]),
    			get_type<Array>(alist_[2]));
    }

    result_type operator()(ArrayFnVariadic fn) const {
    	//std::cout << "array variadic fn" << "\n";
    	std::vector<Array> alist;
    	for(ParserResult item : alist_)
    		alist.push_back(get_type<Array>(item));
    	return fn(alist);
    }

    result_type operator()(IntListFn fn) const {
    	std::vector<int> alist;
    	//std::cout << "int list fn" << "\n";
    	for(ParserResult item : alist_)
    		alist.push_back(get_type<int>(item));
    	return fn(alist);
    }

    result_type operator()(SubscribeFn fn) const {
    	BP_ASSERT(alist_.size() == 2);
    	//std::cout << "subscribe fn" << "\n";
    	return fn(
    			get_type<Array>(alist_[0]),
				get_type<RangeList>(alist_[1]));
    }

    result_type operator()(RangeListFn fn) const {
    	// convert ResultList to RangeList
    	RangeList r_list;
    	//std::cout << "range list fn" << "\n";
    	for(auto item : alist_) {
    		r_list.push_back(get_type<Range>(item));
    	}

    	return fn(r_list);
    }

    result_type operator()(IndexFn fn) const {
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

	result_type operator()(IndexList index_range) const {
		MultiIdxRange res(range_);
		res.sub_range(axis_, index_range);
		return res;
	}

};

struct NamedArrayFn {
	std::string repr;
	ArrayFn fn;
};

inline Range create_index_array(const IndexList& l) {
	return l;
}

inline Range create_slice(const IndexList& l) {
	Slice s;
	std::copy_n(std::make_move_iterator(l.begin()), 3, s.begin());
	return s;
}

inline Range create_index(int idx) {
	return idx;
}

inline Array subscribe(const Array &a, const RangeList &slice_list) {
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

inline RangeList range_list(const RangeList &slice_list) {
	return slice_list;
}

inline Array close_chain(const ComparisonPair & x) {
	return x.first;
}


} // namespace bparser
#endif /* INCLUDE_ARRAY_AST_INTERFACE_HH_ */
