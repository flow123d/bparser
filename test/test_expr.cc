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
	using bparser::expr;

	Array s_const = Const(3.14);
	Array v_const = Const({1, 2, 3});
	Array t_const = Const({{1, 2, 3}, {2, 4, 5}, {3, 5, 6}});

	Array sa = Shape(1);
	Array vb = Shape(3);
	Array tc = Shape(3, 3);

	// scalar expressions
	Array res1 = s_const + sa;
	Array res2 = s_const * sa + sa;
	Array res3 = sa * sa + s_const * sa;

	// elementwise expressions - direct
	Array res10 = vb + vb;
	Array res11 = exp(vb);
	Array res20 = tb + tb;
	Array res21 = exp(tb);

	// broadcasting, slicing
	Array res30 = vb + sa(None());
	Array res31 = tc + vb(None(), Slice());
	Array res32 = vb(Slice(1, None())) +  vb(Slice(0, -1)) + vb(Slice(0, None(), 2));

	// mat mult
	Array res40 = vb.t() ^ vb;
	Array res41 = vb.t() ^ tc ^ vb;
}
