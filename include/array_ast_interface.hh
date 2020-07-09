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
#include "array.hh"

namespace bparser {


typedef std::vector<Array> ArrayList;
typedef std::vector<int> IndexList;
typedef std::array<int, 3> Slice;
typedef boost::variant<IndexList, Slice> Range;
typedef std::vector<Range> RangeList;

typedef boost::variant<
		ArrayList,
		IndexList,
		Slice,
		RangeList
		> ResultList;


struct append_visitor : boost::static_visitor<ResultList> {
	ResultList head_;

	explicit append_visitor(ResultList head)
	: head_(head)
	{}


	Slice operator()(const Slice& alist) const {
		BP_ASSERT(alist.size() == 1);
		ResultList res(head_);
		if (Slice *r = boost::get<Slice>(&res)) {
			return *r;
		} else {
			BP_ASSERT(false);
		}
	}


	template<class T>
	ResultList operator()(const T& alist) const {
		BP_ASSERT(alist.size() == 1);
		ResultList res(head_);
		if (T *r = boost::get<T>(&res)) {
			r->push_back(alist[0]);
			return res;
		} else {
			BP_ASSERT(false);
		}
	}
};

struct size_visitor : boost::static_visitor<uint> {
	template<class T>
	uint operator()(const T& alist) const {
		return alist.size();
	}

};

inline uint result_size(ResultList r) {
	return boost::apply_visitor(size_visitor(), r);
}

typedef Array (*ArrayFnUnary)(const Array &);
typedef Array (*ArrayFnBinary)(const Array &, const Array &);
typedef Array (*ArrayFnTernary)(const Array &, const Array &, const Array &);
typedef Array (*ArrayFnVariadic)(const ArrayList &);
typedef Range (*IntListFn)(const IndexList &);
typedef RangeList (*RangeListFn)(const RangeList &);


typedef boost::variant <
	ArrayFnUnary,
	ArrayFnBinary,
	ArrayFnTernary,
	ArrayFnVariadic,
	IntListFn,  //IndexArray and Slice
	RangeListFn  // slicing
	> ArrayFn;

struct call_visitor {
	typedef ResultList result_type;
	const result_type &alist_;
    explicit call_visitor(const result_type &alist)
	: alist_(alist)
	{}

    template <class List>
    static List get_list(result_type r_list) {
    	if (List* l = boost::get<List>(&r_list)) {
    		return *l;
    	} else {
    		BP_ASSERT(false);
    	}
    }

    result_type operator()(ArrayFnUnary fn) const {
    	BP_ASSERT(result_size(alist_) == 1);
    	ArrayList al = get_list<ArrayList>(alist_);
    	Array res = fn(al[0]);
    	return ArrayList({res});
    }

    result_type operator()(ArrayFnBinary fn) const {
    	BP_ASSERT(result_size(alist_) == 2);
    	ArrayList al = get_list<ArrayList>(alist_);
    	Array res = fn(al[0], al[1]);
    	return ArrayList({res});
    }

    result_type operator()(ArrayFnTernary fn) const {
    	BP_ASSERT(result_size(alist_) == 3);
    	ArrayList al = get_list<ArrayList>(alist_);
    	Array res = fn(al[0], al[1], al[2]);
    	return ArrayList({res});
    }

    result_type operator()(ArrayFnVariadic fn) const {
    	return ArrayList({fn(get_list<ArrayList>(alist_))});
    }

    result_type operator()(IntListFn fn) const {
   		return RangeList({fn(get_list<IndexList>(alist_))});
    }

    result_type operator()(RangeListFn fn) const {
   		return fn(get_list<RangeList>(alist_));
    }
};



// boost::apply_visitor(list_init_visitor
//struct list_init_visitor  : public boost::static_visitor<> {
//	template <class T>
//	T operator()(T x) {
//
//	}
//};

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


}

#endif /* INCLUDE_ARRAY_AST_INTERFACE_HH_ */
