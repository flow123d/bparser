/*
 * nitpick_def.c
 *
 *  Created on: Dec 7, 2025
 *      Author: LV
 */

/*
* Prepare an environment for both DAG generation and then subsequent running
*/

#ifndef NITPICK_INCLUDE_ONLY
#include "nitpick_include.hh"
#endif //NITPICK_INCLUDE_ONLY

using namespace bparser;

#ifndef NITPICK_INCLUDE_ONLY
int main() {
#endif //NITPICK_INCLUDE_ONLY

	// Define own value vectors, preferably aligned.
	constexpr uint vec_size = 8;

	double v1[vec_size * 3];
	for (uint i = 0; i < vec_size * 3; ++i) {
		v1[i] = i;
	}
	double v2[vec_size * 3];
	for (uint i = 0; i < vec_size * 3; ++i) {
		v2[i] = 2;
	}
	double vres[vec_size * 3];

	// Create parser, give the size of the value spaces.
	// That is maximal allocated space. Actual values and 
	// active subset can be changed between evaluations of the expression.
	Parser p(vec_size);
	// parse an expression.
	p.parse("1 * v1 + cs1 * v2");

	// "cs1" constant with shape {}, i.e. scalar and values {2}.
	p.set_constant("cs1", {}, { 2 });
	// "cv1" vector constant with shape {3}
	p.set_constant("cv1", { 3 }, { 1, 2, 3 });
	// "v1" variable with shape {3}; v1 is pointer to the value space
	p.set_variable("v1", { 3 }, v1);
	p.set_variable("v2", { 3 }, v2);
	// Set the result variable (the last line of the expression)
	p.set_variable("_result_", { 3 }, vres);

#ifndef NITPICK_INCLUDE_ONLY
}
#endif //NITPICK_INCLUDE_ONLY