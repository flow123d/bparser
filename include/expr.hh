/*
 * expr.hh
 *
 *  Created on: Dec 26, 2019
 *      Author: jb
 */

#ifndef INCLUDE_EXPR_HH_
#define INCLUDE_EXPR_HH_

#include <memory>
#include <vector>
#include "config.hh"
#include "scalar_expr.hh"

namespace bparser {
namespace expr {
using namespace ::bparser::details;


/**
 * Array nodes represents array operations (with array results) visible to user.
 * These are stored as Array of ScalarNodes, there are no derived classes to ArrayNode which simplifies implementation of operators.
 * We implement overloaded operators and other syntactic sugar in order to create the tree from the C++ code as well as from the parser.
 * This also makes operations from parser more readable.
 *
 * Operations:
 *  - elementwise: functions, operators, use common algorithm for expansion to the scalar nodes using related scalar nodes
 *  - stack - make a higher order from list of lower order
 *    stick - append new "column"
 *  - explicit broadcasting -
 */
struct Array {
	typedef ScalarNode * ScalarNodePtr;

//	typedef std::vector<int> VecInt;
//	typedef std::vector<uint> VecUint;
//
//	/**
//	 * Base set of indices.
//	 */
//	struct IndexSubset {
//		virtual std::vector<uint> idx_list(uint size) const = 0;
//	};
//
//	const int none_int = std::numeric_limits<int>::max;
//
//	// Convertible to none index or none index set
//	struct None {
//		operator int() const {
//			return none_int;
//		}
//
//		operator IndexSubset() const {
//			return IndexList({});
//		}
//	};
//	const None none;
//
//
//	/**
//	 * Index set give by its list.
//	 */
//	struct IndexList : IndexSubset {
//		IndexList(const VecUint &indices)
//		: indices_(indices)
//		{}
//
//		IndexList(std::initializer_list<uint> indices)
//		: indices_(indices)
//		{}
//
//		virtual VecUint idx_list(uint size) const {
//			return indices_;
//		}
//	private:
//		VecUint indices_;
//	};
//
//	/**
//	 * Convertible to the index subset provided size
//	 */
//	struct Slice : IndexSubset {
//		Slice(int begin = 0, int end = none, int step = 1)
//		: begin_(begin),
//		  end_(end),
//		  step_(step)
//		{}
//
//		VecUint idx_list(uint size) const override {
//			VecUint indices;
//			if (begin_ < 0) begin_ += size;
//			if (end_ == int(none)) end_ = size;
//			if (end_ < 0) end_ += size;
//			for(int i = begin_; (end_ - begin_) * step_ > 0; i+=step_)
//				indices.push_back(i);
//			return std::move(indices);
//		}
//
//		int begin_, end_, step_;
//	};
//
//	struct MultiIdx {
//		MultiIdx(const VecUint &indices)
//		: indices_(indices)
//		{}
//
//		bool operator!=(const MultiIdx &other) {
//			return indices_ != other.indices_;
//		}
//
//		bool inc(const VecUint &shape) {
//			ASSERT(shape.size() == indices_.size());
//			for(uint i = 0; i < shape.size(); ++i) {
//				indices_[i] += 1;
//				if (indices_[i] == shape[i]) {
//					indices_[i] = 0;
//					continue;
//				} else {
//					return true;
//				}
//			}
//			return false;
//		}
//
//		VecUint indices_;
//	};
//
//	struct MultiIdxRange {
//		MultiIdxRange(const std::vector<IndexSubset> &range_spec) {
//			for(IndexSubset subset : range_spec) ranges_.push_back(subset.idx_list());
//		}
//
//		MultiIdx begin();
//		MultiIdx end();
//		std::vector<std::vector<uint>> ranges_;
//	};
//


	/**
	 * ?? Do we need default constructor?
	 */
	Array()
	{}



	// Const scalar node.
	Array(double x)
	: elements_({Scalarode}){

	}

	// Const 1D array node.
	Array(std::initializer_list<double> list);

	// Const 2D array node.
	Array(std::initializer_list<std::initializer_list<double>> list);

	static Array stick(Array column);

	static Array broadcast();
	/**
	 * Need indexing syntax, preferably close to the Python syntax selected for the parser:
	 *
	 * Python				C++ interface
	 * a[0, 1]      		a(0, 1)
	 * a[[0,1,3], 3]		a({0,1,3}, 3)
	 * a[:, 3]				a(Slice(), 3)
	 * a[1:, 3]				a(Slice(1), 3)
	 * a[:-1, 3]		    a(Slice(0, -1), 3)
	 * a[::2, 3]		    a(Slice(0, None(), 2), 3)
	 *
	 * a[None, :]			a(None(), Slice())
	 *
	 *
	 * That is slice() is special, but then we haove only slice(start=0,
	 *
	 * TODO:
	 * 1. variadic arguments for the operator() have to be converted to a vector of IndexSusets
	 * 2. same we get from the parser
	 * 3. Index Subsets (with known shape) converted to the vector of index sets, the we use the MultiIdx to iterate over
	 *  own and other Array.
	 *
	 *  Not clear how to apply operations.
	 */

	 void set_subarray_(MultiIdx idx_out, MultiIdx idx_in, Array other, in_subset, ...) {

	 }

private:
	ScalarExpression se;
	std::vector<scalar_node_ptr> elements_;
	std::vector<uint> shape_;

};

} // expr
} // bparser


#endif /* INCLUDE_EXPR_HH_ */
