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
#include <algorithm>
#include <cmath>
#include <boost/math/constants/constants.hpp>
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
	typedef std::vector<uint> Shape;

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

	struct MultiIdxRange {

		MultiIdxRange(Shape shape)
		: full_shape_(shape) {
			// full range
			for(uint axis=0; axis < full_shape_.size(); ++axis) {
				std::vector<uint> full(full_shape_[axis]);
				for(uint i=0; i < full_shape_[axis]; ++i) full[i] = i;
				ranges_.push_back(full);
			}
		}

//		void subset(const std::vector<IndexSubset> &range_spec) {
//			ASSERT(range_spec.size() == full_shape_.size());
//			for(IndexSubset subset : range_spec)
//				ranges_.push_back(subset.idx_list(full_shape_));
//		}


		/**
		 * Broadcast the full range to the given shape.
		 * 1. pad full shape from left by ones
		 * 2. extend dimensions equal to one to the given shape
		 * 3. other dimensions must match given shape
		 * TODO: make all methods either pure or inplace
		 */
		MultiIdxRange broadcast(Shape other) {
			int n_pad = other.size() - full_shape_.size();
			if ( n_pad < 0) {
				std::ostringstream ss;
				ss << "Broadcast from longer shape " << full_shape_.size() <<
						" to shorter shape " << other.size() << ".\n";
				Throw(ss.str());
			}
			Shape res_shape(n_pad, 1);
			res_shape.insert(res_shape.end(), full_shape_.begin(), full_shape_.end());

			MultiIdxRange result(other);
			for(uint ax=0; ax < res_shape.size(); ++ax) {
				if (res_shape[ax] == 1)
					result.ranges_[ax] = std::vector<uint>(other[ax], 0);
				else if (res_shape[ax] != other[ax]) {
					std::ostringstream ss;
					ss << "Broadcast from " << res_shape[ax] << " to "
					   << other[ax] << " in axis " << ax;
					Throw(ss.str());
				}

			}
			return result;
		}

		static Shape broadcast_common_shape(Shape a, Shape b) {
			uint res_size = std::max(a.size(), b.size());
			Shape res(res_size);
			a.insert(a.begin(), res_size - a.size(), 1);
			b.insert(b.begin(), res_size - b.size(), 1);
			for(uint i=0; i<res_size; ++i) {
				if (a[i]==1) a[i] = b[i];
				if (b[i]==1) b[i] = a[i];
				if (a[i] == b[i]) {
					res[i] =  a[i];
				} else {
					std::ostringstream ss;
					ss << "Common broadcast between " << a[i] << " and "
					   << b[i] << " in axis " << i;
					Throw(ss.str());
				}

			}
			return res;
		}

		// Insert 'axis' with 'dimension' into range
		// E.g. for axis=1, dimension=2:
		// range [[1,2],[2,5,6]] -> [[1,2],[0,1],[2,5,6]]
		void insert_axis(uint axis, uint dimension=1) {
			full_shape_.insert(full_shape_.begin() + axis, dimension);
			ranges_.insert(ranges_.begin() + axis, std::vector<uint>());
			for(uint i=0; i<dimension; ++i)
				ranges_[axis].push_back(i);
		}

		// Shift the range in given 'axis' by given 'shift'.
		void shift_axis(uint axis, uint shift) {
			for(uint &el : ranges_[axis]) el +=shift;
		}

		/**
		 * repeat last element in 'axis' range
		 */
		void pad(uint n_items, uint axis) {
			ranges_[axis].insert(ranges_[axis].end(), n_items, ranges_[axis].back());
		}


		/**
		 * Range resulting shape.
		 */
		Shape sub_shape() const {
			Shape shape(ranges_.size());
			for(uint axis=0; axis < shape.size(); ++axis)
				shape[axis] = ranges_[axis].size();
			return shape;
		}

		std::vector<std::vector<uint>> ranges_;
		Shape full_shape_;
	};

	struct MultiIdx {
		typedef std::vector<uint> VecUint;

		MultiIdx(MultiIdxRange range)
		: range_(range), indices_(range.full_shape_.size(), 0)
		{}


		bool inc() {
			for(uint axis = indices_.size(); axis > 0 ; )
			{
				--axis;
				indices_[axis] += 1;
				if (indices_[axis] == range_.ranges_[axis].size()) {
					indices_[axis] = 0;
					continue;
				} else {
					return true;
				}
			}
			return false;
		}

		// Linear index into Array::elements_.
		// Last index runs fastest.
		uint linear_idx() {
			uint lin_idx = 0;
			for(uint axis = 0; axis < indices_.size(); ++axis) {
				lin_idx *= range_.full_shape_[axis];
				lin_idx += range_.ranges_[axis][indices_[axis]];
			}
			return lin_idx;
		}

		uint linear_subidx() {
			uint lin_idx = 0;
			for(uint axis = 0; axis < indices_.size(); ++axis) {
				lin_idx *= range_.ranges_[axis].size();
				lin_idx += indices_[axis];
			}
			return lin_idx;
		}



		// List of indices for every axis.
		MultiIdxRange range_;
		// Indirect multiindex index into range_.
		// Last index runs fastest.
		VecUint indices_;
	};



	MultiIdxRange range() const {
		return MultiIdxRange(shape_);
	}

	ScalarNodePtr &operator[](MultiIdx idx) {
		return elements_[idx.linear_idx()];
	}

	ScalarNodePtr operator[](MultiIdx idx) const {
		return elements_[idx.linear_idx()];
	}

	//Array operator()()

	static Array deg_to_rad_factor() {
		static ConstantNode f( boost::math::constants::pi<double>() / 180 );
		Array a;
		a.elements_[0] = &f;
		return a;
	}

	static Array rad_to_deg_factor() {
		static ConstantNode f(  180 / boost::math::constants::pi<double>() );
		Array a;
		a.elements_[0] = &f;
		return a;
	}

	template <class T>
	static Array unary_op(const Array &a) {
		Array result(a.shape_);
		MultiIdx idx(a.range());
		for(;;) {
			result[idx] = ScalarNode::create<T>(a[idx]);
			if (! idx.inc()) break;
		}
		return result;
	}


	template <class T>
	static Array binary_op(const Array &a, const Array &b) {

		Shape res_shape = MultiIdxRange::broadcast_common_shape(a.shape_, b.shape_);
		MultiIdxRange a_range = a.range().broadcast(res_shape);
		MultiIdxRange b_range = b.range().broadcast(res_shape);

		MultiIdx a_idx(a_range);
		MultiIdx b_idx(b_range);
		Array result(res_shape);
		for(;;) {
			ASSERT(a_idx.linear_subidx() == b_idx.linear_subidx());
			result.elements_[a_idx.linear_subidx()] =
					ScalarNode::create<T>(
							a.elements_[a_idx.linear_idx()],
							b.elements_[b_idx.linear_idx()]);
			if (!a_idx.inc() || !b_idx.inc()) break;
		}
		return result;
	}



	/**
	 * ?? Do we need default constructor?
	 */
	Array()
	{}

	Array(const Array &other)
	: Array(other, other.range())
	{}

	Array operator=(const Array & other)
	{
		elements_ = other.elements_;
		shape_ = other.shape_;
		return *this;
	}


	Array(const Shape &shape)
	: shape_(shape)
	{
		uint shape_prod = 1;
		for(uint dim : shape) shape_prod *= dim;
		elements_.resize(shape_prod);
	}

	/**
	 * Subshape constructor.
	 */
	Array(const Array &other, const MultiIdxRange &range)
	: Array(range.sub_shape())
	{
		from_subset(other, range);
	}

	/// Create a result array from *this using storage of 'variable'.
	Array make_result(const Array &variable) {
		Array res(shape_);
		for(uint i=0; i<elements_.size(); ++i)
			res.elements_[i] = ScalarNode::create_result(elements_[i], variable.elements_[i]->values_);
		return res;
	}

	const std::vector<ScalarNodePtr> &elements() const {
		return elements_;
	}


	Shape minimal_shape(Shape other) {
		Shape result;
		for(uint e : other)
			if (e > 1) result.push_back(e);
		return result;
	}

	/**
	 * Set subset of *this given by the range to values given by 'other'.
	 * The shape of range have length same as 'shape_' but non 1 dimensions have to match
	 * shape of other.
	 */
	void set_subset(const MultiIdxRange &range, const Array &other) {
		Shape sub_shp = range.sub_shape();
		Shape min_shape_a = minimal_shape(sub_shp);
		Shape min_shape_b = minimal_shape(other.shape_);
		ASSERT(min_shape_a == min_shape_b);
		MultiIdx idx(range);
		for(;;) {
			elements_[idx.linear_idx()] = other.elements_[idx.linear_subidx()];
			if (! idx.inc()) break;
		}
	}

	/**
	 * Fill all elements of *this to the subset of 'other' given by the 'range'.
	 */
	void from_subset(const Array &other, const MultiIdxRange &range) {
		MultiIdx idx(range);
		for(;;) {
			elements_[idx.linear_subidx()] = other.elements_[idx.linear_idx()];
			if (! idx.inc()) break;
		}
	}




	void _set_shape(Shape shape) {
		shape_ = shape;
		uint shape_prod = 1;
		for(uint dim : shape) shape_prod *= dim;
		elements_.resize(shape_prod);
	}



	// Const scalar node.
	static Array constant(const std::vector<double> &values, Shape shape = {}) {
		Array res(shape);
		for(uint i_el=0; i_el < res.elements_.size(); ++i_el)
			res.elements_[i_el] = ScalarNode::create_const(values[i_el]);
		return res;
	}


//	// Const 1D array node.
//	// TODO: do better possibly without transform it is unable to resolve overloaded methods
//	static Array constant1(const std::vector<double> &list)
//	{
//		std::vector<Array> ar_list;
//		std::transform(list.begin(), list.end(), std::back_inserter(ar_list), Array::constant);
//		return Array::stack(ar_list);
//	}
//
//	// Const 2D array node.
//	static Array constant2(const std::vector<std::vector<double>> &list)
//	{
//		std::vector<Array> ar_list;
//		std::transform(list.begin(), list.end(), std::back_inserter(ar_list), Array::constant1);
//		return Array::stack(ar_list);
//	}


	// Vector value array with given shape
	static Array value(double *v, uint array_max_size, Shape shape = {})
	{
		Array res(shape);
		for(uint i_el=0; i_el < res.elements_.size(); ++i_el) {
			res.elements_[i_el] = ScalarNode::create_value(v);
			v += array_max_size;
		}
		return res;

//		uint i_el=0;
//		Shape idx(shape.size(), 0);
//		int ax = idx.size() - 1;
//		for(;;) {
//			res.elements_[i_el++] = ScalarNode::create_value(v);
//			v += array_max_size;
//			for(;;) {
//				if (ax < 0) goto exit;
//				idx[ax]++;
//				if (idx[ax] < shape[ax]) {
//					ax = 0;
//					break;
//				} else {
//					idx[ax] = 0;
//					ax--;
//				}
//			}
//		}
//		exit:
//		return res;
	}



	/**
	 * 1. broadcast slice to have same shape in other axes then 'axis'
	 * 2.
	 */
	Array append(const Array &slice, uint axis=0)
	{
		return Array::stick({*this, slice}, axis);
	}


	static Array stack(const std::vector<Array> &list, uint axis=0) {
		if (list.size() == 0) throw;
		if (list.size() == 1) return list[0];
		// same shape
		for(auto ar : list) ar.shape_ == list[0].shape_;

		// new shape
		Shape res_shape = list[0].shape_;
		if (axis < 0) axis = res_shape.size() + axis;
		ASSERT(axis >= 0);
		res_shape.insert(res_shape.begin() + axis, list.size());

		// 3. broadcast arrays and concatenate
		Array result = Array(res_shape);

		for(uint ia=0; ia < list.size(); ++ia) {
			MultiIdxRange bcast_range = list[0].range();
			bcast_range.insert_axis(axis);
			bcast_range.shift_axis(axis, ia);
			result.set_subset(bcast_range, list[ia]);
		}
		return result;

	}


	static Array stick(const std::vector<Array> &list, uint axis=0) {
		if (list.size() == 0) throw;
		if (list.size() == 1) return list[0];

		// TODO: iteratively prepare multi index ranges (can be modified in place)
		// then apply at once
		// Ranges - here we need to construct subranges for the new array
		// TODO:

		// determine resulting shape
		// 1. max shape length
		size_t max_len=0;
		for(auto ar : list) max_len = std::max(max_len, ar.shape_.size());
		// 2. max shape
		Shape res_shape(max_len, 1);
		std::vector<MultiIdxRange> list_range;
		uint axis_dimension = 0;
		for(auto ar : list) {
		    MultiIdxRange r = ar.range();
		    uint sub_size = r.ranges_[axis].size();
		    if ( sub_size < max_len -1) throw;
		    if ( sub_size < max_len) {
		    	r.insert_axis(axis);
		    }
			for(uint ax=0; ax < res_shape.size(); ++ax)
				res_shape[ax] = std::max(res_shape[ax], (uint)(r.ranges_[ax].size()));
			axis_dimension += r.ranges_[axis].size();
			list_range.push_back(r);
		}

		// final shape, sum along axis
		// first range, pad to final shape
		res_shape[axis] = axis_dimension;
		// 3. broadcast arrays and concatenate
		Array result = Array(res_shape);
		for(uint ia=1; ia < list.size(); ++ia) {
			Shape s = list_range[0].sub_shape();
			s[axis] = list[0].shape_[axis];
			MultiIdxRange r = list_range[0].broadcast(s);
			r.shift_axis(axis, 0);
			result.set_subset(r, list[0]);
		}
		return result;
	}



	/**
	 * - if this->shape_[i] == 1 the array is repeated in that axis according to other[i]
	 * - this->shape_ is padded from the left by ones
	 * - final array must have same shape as other
	 */
	Array broadcast(Shape other) {
		MultiIdxRange bcast_range = range().broadcast(other);
		return Array(*this, bcast_range);
	}

	/**
	 * Need indexing syntax, preferably close to the Python syntax selected for the parser:
	 *
	 * Python				C++ interface
	 * a[0, 1]      		a(0, 1)
	 * a[[0,1,3], 3]		a({0,1,3}, 3)
	 * a[:, 3]				a(Slice(), 3)
	 * a[1:, 3]				a(Slice(1), 3)
	 * a[:-1, 3]		    a(Slice(0, -1), 3)
	 * a[::2, 3]		    a(Slice(0, None, 2), 3)
	 *
	 * a[None, :]			a(None, Slice())
	 *
	 *
	 *
	 * TODO:
	 * 1. variadic arguments for the operator() have to be converted to a vector of IndexSubsets
	 * 2. same we get from the parser
	 * 3. Index Subsets (with known sub_shape) converted to the vector of index sets, the we use the MultiIdx to iterate over
	 *  own and other Array.
	 *
	 *  Not clear how to apply operations.
	 */

//	 void set_subarray_(MultiIdx idx_out, MultiIdx idx_in, Array other, in_subset, ...) {
//
//	 }

private:
	std::vector<ScalarNodePtr> elements_;
	Shape shape_;

};

template <class Fn>
Array func(const Array &x) {
	return Array::unary_op<Fn>(x);
}

/**
 * Special function
 */
Array deg_fn(const Array & rad) {
	return Array::binary_op<_mul_>(rad, Array::deg_to_rad_factor());
}

Array rad_fn(const Array & deg) {
	return Array::binary_op<_mul_>(deg, Array::rad_to_deg_factor());
}

expr::Array gt_op(const expr::Array & a, const expr::Array & b) {
	return Array::binary_op<_lt_>(b, a);
}

expr::Array ge_op(const expr::Array & a, const expr::Array & b) {
	return Array::binary_op<_le_>(b, a);
}

expr::Array unary_plus(const expr::Array & a) {
	return a;
}


Array operator+(const Array &a,  const Array &b) {
	return Array::binary_op<_add_>(a, b);
}

Array operator*(const Array &a,  const Array &b) {
	return Array::binary_op<_mul_>(a, b);
}

} // expr
} // bparser


#endif /* INCLUDE_EXPR_HH_ */
