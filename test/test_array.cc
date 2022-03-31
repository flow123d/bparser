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
	EXPECT_CONTINUE(TEST_EQ(shape_size(r.target_shape()), sub_linear.size()));
	for(uint n=0;n < sub_linear.size();n++) {
		std::cout << print_vector(idx.src_indices_) << "\n";
		std::cout << n << ", " << idx.idx_src() << ", " << sub_linear[n] <<  ", " << idx.idx_trg() << "\n";
		EXPECT_CONTINUE(TEST_EQ(n, idx.idx_trg()));
		EXPECT_CONTINUE(TEST_EQ(sub_linear[n], idx.idx_src()));
		if (! idx.inc_trg()) break;
	}
	EXPECT_CONTINUE(TEST_EQ(idx.idx_trg(), uint(0))); // all values compared
	return true;
}

bool shape_eq(Shape a, Shape b) {
	bool match = vec_eq(a, b);
	if (!match) std::cout << "Different shapes!\n\n";
	return match;}

void test_constructors() {
	Array a;
	EXPECT(a.is_none());
	EXPECT(TEST_EQ(a.range().source_shape_.size(), uint(0)));
	Array b;
	b=a;
	EXPECT(b.is_none());
	Array c(a);
	EXPECT(c.is_none());
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
	{
		auto r = MultiIdxRange(Shape({2,1})).full();
		EXPECT(match_linear(r, {0,1}));

		// broadcasting and MultiIdx test
		EXPECT(shape_eq(MultiIdxRange::broadcast_common_shape({5,1,3}, {2,1}), {5, 2, 3}) );
		MultiIdxRange bcast_r = r.broadcast(Shape({2,2,3}));
		EXPECT(shape_eq(bcast_r.source_shape_, {2, 1}));
		EXPECT(shape_eq(bcast_r.target_shape(), {2, 2, 3}));
		EXPECT(shape_eq(bcast_r.ranges_[0], {0,0}));
		EXPECT(shape_eq(bcast_r.ranges_[1], {0,1}));
		EXPECT(shape_eq(bcast_r.ranges_[2], {0, 0, 0}));
		EXPECT(shape_eq(bcast_r.target_transpose_, {0, 1, 2}));
		EXPECT(match_linear(bcast_r, {0,0,0,1,1,1,0,0,0,1,1,1}));
		ASSERT_THROW(r.broadcast(Shape({4,3,3})),
							"Invalid broadcast from 2 to 3");
	}
	// test MIR modifications
	{
		auto r = MultiIdxRange(Shape({2,1})).full();
		//insert_axis
		r.insert_axis(0, 0, 3);
		EXPECT(shape_eq(r.target_shape(), {3, 2, 1}));
		EXPECT(shape_eq(r.source_shape_, {3, 2, 1}));
	}

	// subscriptions
	{
		auto r = MultiIdxRange(Shape({2,3,2}));
		//insert_range

		r.sub_range({1,0}); // axis 0
		ASSERT_THROW(r.sub_range({3}), "out of range");
		r.sub_slice({none_int, 2, none_int}); // axis 1
		ASSERT_THROW(r.sub_slice({none_int, none_int, 0}),
				"Slice step cannot be zero.");
		r.sub_none(); // axis 2

		r.sub_index(1); // dest axis 3 / src axis 2
		EXPECT(shape_eq(r.target_shape(), {2,2,1}))
		EXPECT(match_linear(r, {7,9, 1,3}));

	}
	{
		auto r = MultiIdxRange(Shape({2,3,2}));
		//insert_range

		r.sub_range({1});
		r.sub_slice({2, none_int,-2});
		r.sub_index(0);
		EXPECT(match_linear(r, {10,6}));
		EXPECT(shape_eq(r.target_shape(), {1,2}))
	}
	{
		// swap axis
		auto r = MultiIdxRange(Shape({2,3,2}));
		r.sub_range({1,0}); // axis 0
		r.sub_slice({none_int, 2, none_int}); // axis 1
		r.sub_none(); // dest axis 2
		r.sub_index(1); //  src axis 2

		r.target_transpose(-1, -2);
		EXPECT(shape_eq(r.source_shape_ , {2,3,2}));
		EXPECT(shape_eq(r.target_shape() , {2,1,2}));
		EXPECT(vec_eq(r.ranges_[0],{1,0}));
		EXPECT(vec_eq(r.ranges_[1],{0,1}));
		EXPECT(vec_eq(r.ranges_[2],{1}));
		EXPECT(vec_eq(r.target_transpose_, {0, uint(none_int), 1}));

	}
	{
		// swap axis
		MultiIdxRange r = MultiIdxRange({2,1}).broadcast({3,2,1});
		r.target_transpose(-2, -1);
		EXPECT(shape_eq(r.source_shape_, {2,1}));
		EXPECT(shape_eq(r.target_shape() , {3,1,2}));
		EXPECT(vec_eq(r.ranges_[0],{0,0,0}));
		EXPECT(vec_eq(r.ranges_[1],{0,1}));
		EXPECT(vec_eq(r.ranges_[2],{0}));
		EXPECT(vec_eq(r.target_transpose_, {0,2,1}));
		EXPECT(match_linear(r, {0,1,0,1,0,1}));

	}


}

void test_MultiIdx() {
	// test inc_dest
	{
		auto r = MultiIdxRange(Shape({2,3,2}));
		r.sub_range({1,0}); // axis 0
		r.sub_slice({none_int, 2, none_int}); // axis 1
		r.sub_none(); // axis 2
		r.sub_index(1); // dest axis 3 / src axis 2

		MultiIdx idx(r);
		std::cout << "MultiIdx 0\n";
		EXPECT(vec_eq(idx.indices(), {0,0,0}));
		EXPECT( idx.idx_src() == 7);
		EXPECT( idx.idx_trg() == 0);
		idx.inc_trg();

		std::cout << "MultiIdx 1\n";
		EXPECT(vec_eq(idx.indices(), {0,1,0}));
		EXPECT( idx.idx_src() == 9);
		EXPECT( idx.idx_trg() == 1);
		idx.inc_trg(-1, 2);

		std::cout << "MultiIdx 2\n";
		EXPECT(vec_eq(idx.indices(), {1,1,0}));
		EXPECT( idx.idx_src() == 3);
		EXPECT( idx.idx_trg() == 3);
		idx.inc_trg(-3, -1);

		std::cout << "MultiIdx 4\n";
		EXPECT(vec_eq(idx.indices(), {0,1,0}));
		EXPECT( idx.idx_src() == 9);
		EXPECT( idx.idx_trg() == 1);
		EXPECT( idx.valid());

		std::cout << "MultiIdx 5\n";
		EXPECT( false == idx.inc_trg(1, -2, false));

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

	test_constructors();
	test_absolute_idx();
	test_MultiIdxRange();
	test_MultiIdx();

	test_Array();
}
