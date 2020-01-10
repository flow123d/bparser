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
#include "expr.hh"




int main() {
	using namespace bparser::expr;

	// constant array, array of ConstNode
	Array s_const = Array::constant(3.14);
	Array v_const = Array::constant1(std::vector<double>({1.0, 2, 3}));
	//Array t_const = Array::constant({{1, 2, 3}, {2, 4, 5}, {3, 5, 6}});

	double v[20][3][3];
	// variable arrays, array of ValueNode
	Array sa = Array::value((double *)v, 20); // double* ptr, shape
	Array vb = Array::value((double *)v, 20, {3});
	Array tc = Array::value((double *)v, 20, {3,3});

//	// constant - vector expressions (no difference)
//	Array res1 = s_const + sa;
//	Array res2 = s_const * sa + sa;
//	Array res3 = sa * sa + s_const * sa;
//
//	// elementwise expressions - direct
//	Array res10 = vb + vb;
//	Array res11 = exp(vb);
//	Array res20 = tb + tb;
//	Array res21 = exp(tb);
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
