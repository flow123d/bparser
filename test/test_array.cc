/*
 * test_array.hh
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */


#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>

#include "array.hh"
#include "test_tools.hh"

using namespace bparser;
using namespace bparser::details;


#define EXPECT_CONTINUE(expr) \
		if (!expr) return false

/**
 * multi_idx.linear() == MIR
 */
bool match_linear(MultiIdxRange r, std::vector<uint> sub_linear) {
	MultiIdx idx(r);
	EXPECT_CONTINUE(TEST_EQ(shape_size(r.sub_shape()), sub_linear.size()));
	for(uint n=0;n < sub_linear.size();n++) {
		//std::cout << print_vector(idx.indices_) << "\n";
		//std::cout << n << ", " << idx.linear_subidx() << ", " << sub_linear[n] <<  ", " << idx.linear_idx() << "\n";
		EXPECT_CONTINUE(TEST_EQ(n, idx.dest_idx()));
		EXPECT_CONTINUE(TEST_EQ(sub_linear[n], idx.src_idx()));
		if (! idx.inc()) break;
	}
	EXPECT_CONTINUE(TEST_EQ(idx.dest_idx(), uint(0))); // all values compared
	return true;
}

bool shape_eq(Shape a, Shape b) {

	TEST_EQ(a.size(), b.size());
	uint size = a.size();
	bool match = true;
	for(uint i=0;i<size; i++)
		if (a[i] != b[i]) {
			match = false;
			break;
		}

	if (!match) {
		std::cout << "Different shapes: \n"
				<< print_vector(a) << "\n"
				<< print_vector(b) << "\n";
		return false;
	}
	return true;
}


void test_absolute_idx() {
	EXPECT(TEST_EQ(uint(0), absolute_idx(0, 5)));
	EXPECT(TEST_EQ(uint(4), absolute_idx(4, 5)));
	ASSERT_THROW(absolute_idx(5, 5),
			"out of range");
	EXPECT(TEST_EQ(uint(4), absolute_idx(-1, 5)));
	EXPECT(TEST_EQ(uint(0), absolute_idx(-5, 5)));
	ASSERT_THROW(absolute_idx(-6, 5),
			"out of range");
}

void test_MultiIdxRange() {
	// TODO: test MultiRange, MultiIdx
	{
		auto r = MultiIdxRange(Shape({2,1})).full();
		EXPECT(match_linear(r, {0,1}));

		// broadcasting and MultiIdx test
		EXPECT(shape_eq(MultiIdxRange::broadcast_common_shape({3}, {2,1}), {2, 3}) );
		MultiIdxRange bcast_r = r.broadcast(Shape({2,2,3}));
		EXPECT(bcast_r.full_shape_ == Shape({2, 2, 3}));
		EXPECT(match_linear(bcast_r, {0,0,0,3,3,3,0,0,0,3,3,3}));
		ASSERT_THROW(r.broadcast(Shape({4,3,3})),
							"Broadcast from 2 to 3");
	}
	// test MIR modifications
	{
		auto r = MultiIdxRange(Shape({2,1})).full();
		//insert_axis
		r.insert_axis(0, 0, 3);
		EXPECT(shape_eq(r.sub_shape(), {3, 2, 1}));
		EXPECT(shape_eq(r.full_shape_, {3, 2, 1}));
	}
	{
		auto r = MultiIdxRange(Shape({2,3,2})).full();
		//insert_range

		r.sub_range(0, 0, {1,0});
		EXPECT(match_linear(r, {6,7,8,9,10,11, 0,1,2,3,4,5}));
		r.sub_range(1, 1, {-2, -1});
		EXPECT(match_linear(r, {8,9,10,11, 2,3,4,5}));
		ASSERT_THROW(r.sub_range(1, 1, {3}),
				"out of range");

		r.sub_slice(1, 1, {none_int, 2, none_int});
		EXPECT(match_linear(r, {6,7,8,9, 0,1,2,3}));
		r.sub_slice(1, 1, {2, none_int,-2});
		EXPECT(match_linear(r, {10,11,6,7, 4,5,0,1}));
		ASSERT_THROW(r.sub_slice(1, 1, {none_int, none_int, 0}),
				"Slice step cannot be zero.");

		r.sub_index(2, 2, 1);
		EXPECT(match_linear(r, {11,7, 5,1}));
	}

}

void test_Array() {
	// constant array, array of ConstNode
	Array s_const = Array::constant({3.14});
	Array v_const = Array::constant({1.0, 2, 3}, {3});
	Array t_const = Array::constant({1, 2, 3, 2, 4, 5, 3, 5, 6}, {3, 3});

	double v[20][3][3];
	// variable arrays, array of ValueNode
	Array sa = Array::value((double *)v, 20); // double* ptr, shape
	Array vb = Array::value((double *)v, 20, {3});
	Array tc = Array::value((double *)v, 20, {3,3});

	// TODO: implement and test equality function
	// TODO: implement and test indexing and slicing

	// constant - vector expressions (no difference)
	Array res1 = s_const + sa;
	Array res2 = s_const * sa + sa;
	Array res3 = sa * sa + s_const * sa;

	// elementwise expressions - direct
	Array res10 = vb + vb;
	Array res11 = func<_exp_>(vb);
	Array res20 = t_const + tc;
	Array res21 = func<_exp_>(tc);

	Array res31;
	res31 = res31.append(s_const);
	res31 = res31.append(s_const);
	res31 = res31.append(sa);
	BP_ASSERT(res31.shape() == Shape({3}));

//
//	// broadcasting, slicing
//	Array res30 = vb + sa(None());
//	Array res31 = tc + vb(None(), Slice());
//	Array res32 = vb(Slice(1, None())) +  vb(Slice(0, -1)) + vb(Slice(0, None(), 2));
//
//	// mat mult
//	Array res40 = vb.t() ^ vb;
//	Array res41 = vb.t() ^ tc ^ vb;

}

int main() {
	test_absolute_idx();
	test_MultiIdxRange();
	test_Array();
}
