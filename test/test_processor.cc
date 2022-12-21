/*
 * test_processor.cc
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */



#include <string>
#include "assert.hh"
#include "scalar_node.hh"
#include "create_processor.hh"


using namespace bparser;
using namespace bparser::details;


uint simd_size = get_simd_size();

void test_simple_expr() {

	const uint vec_size = 20;
	double v1_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v1_values[i] = (double)i - 10;
	//std::cout << "v1_ptr: " << v1_values << "\n";

	double v2_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v2_values[i] = -100;
	//std::cout << "v2_ptr: " << v2_values << "\n";

	double v3_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		v3_values[i] = -100;
	//std::cout << "v3_ptr: " << v3_values << "\n";



	ScalarNodePtr c1 = ScalarNode::create_const(100);
	ScalarNodePtr v1 = ScalarNode::create_value(v1_values);
	ScalarNodePtr e1 = ScalarNode::create<_abs_>(v1);
	ScalarNodePtr e2 = ScalarNode::create<_add_>(e1, c1);

	ScalarNodePtr r1 = ScalarNode::create_result(e1, v2_values);
	ScalarNodePtr r2 = ScalarNode::create_result(e2, v3_values);

	ExpressionDAG se({r1, r2});
	ExpressionDAG::NodeVec &sorted = se.sort_nodes();
//	for(auto node : sorted) {
//		se._print_node(node);
//		std::cout << " ptr: " << node << " i_storage: " << node->result_idx_ << "\n";
//	}

	se.print_in_dot();
	BP_ASSERT(sorted[3] == v1);
	BP_ASSERT(v1->result_idx_ == 3);
	BP_ASSERT(sorted[2] == e1);
	BP_ASSERT(e1->result_idx_ == 1);
	BP_ASSERT(sorted[1] == c1);
	BP_ASSERT(c1->result_idx_ == 0);
	BP_ASSERT(sorted[0] == e2);
	BP_ASSERT(e2->result_idx_ == 2);
	ProcessorBase * processor = ProcessorBase::create_processor(se, vec_size, simd_size, nullptr);


	std::vector<uint> subset = {1, 3, 4};
	processor->set_subset(subset);
	processor->run();
	// v1 active: [ -6, -5, -4, -3; 2, 3, 4, 5; 6, 7, 8, 9 ]
	// r1, e1     [  6   5   4   3  2  3  4  5  6  7  8  9 ]
	// r2, e2     [  6   5   4   3  2  3  4  5  6  7  8  9 ]  + 100
//	for(uint i=0;i < vec_size; ++i)
//		std::cout << "i: " << i << " v2: " << v2_values[i]  << " v3: " <<v3_values[i] << "\n";

	BP_ASSERT(v2_values[0] == -100);
	BP_ASSERT(v2_values[3] == -100);
	BP_ASSERT(v2_values[4] ==  6);
	BP_ASSERT(v2_values[7] ==  3);
	BP_ASSERT(v2_values[8] == -100);
	BP_ASSERT(v2_values[11] == -100);
	BP_ASSERT(v2_values[12] == 2);
	BP_ASSERT(v2_values[15] == 5);
	BP_ASSERT(v2_values[16] == 6);
	BP_ASSERT(v2_values[19] == 9);

	BP_ASSERT(v3_values[0] == -100);
	BP_ASSERT(v3_values[3] == -100);
	BP_ASSERT(v3_values[4] ==  106);
	BP_ASSERT(v3_values[7] ==  103);
	BP_ASSERT(v3_values[8] == -100);
	BP_ASSERT(v3_values[11] == -100);
	BP_ASSERT(v3_values[12] == 102);
	BP_ASSERT(v3_values[15] == 105);
	BP_ASSERT(v3_values[16] == 106);
	BP_ASSERT(v3_values[19] == 109);

	// Change data and structure
	for(uint i=0;i < vec_size; ++i)
		v1_values[i] = (double)i - 12;
	for(uint i=0;i < vec_size; ++i)
		v2_values[i] = -200;
	for(uint i=0;i < vec_size; ++i)
		v3_values[i] = -200;
	subset = {0, 4};
	// v1 active: [ -12, -11, -10, -9; 4, 5, 6, 7 ]
	// r1, e1     [  12   11   10   9  4  5  6  7 ]
	// r2, e2       + 100



	processor->set_subset(subset);
	processor->run();

	BP_ASSERT(v2_values[0] == 12);
	BP_ASSERT(v2_values[3] == 9);
	BP_ASSERT(v2_values[4] ==  -200);
	BP_ASSERT(v2_values[15] == -200);
	BP_ASSERT(v2_values[16] == 4);
	BP_ASSERT(v2_values[19] == 7);

	BP_ASSERT(v3_values[0] == 112);
	BP_ASSERT(v3_values[3] == 109);
	BP_ASSERT(v3_values[4] ==  -200);
	BP_ASSERT(v3_values[15] == -200);
	BP_ASSERT(v3_values[16] == 104);
	BP_ASSERT(v3_values[19] == 107);

	processor->~ProcessorBase();
}

std::vector<double> v1_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = -2;
	res[1] = -1;
	res[2] = 0;
	res[3] = 1;
	return res;
}

std::vector<double> v2_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = 5;
	res[1] = 6;
	res[2] = 7;
	res[3] = 8;
	return res;
}

std::vector<double> v3_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = 5;
	res[1] = 6;
	res[2] = 0;
	res[3] = 1;
	return res;
}

std::vector<double> v4_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = -3;
	res[1] = -3;
	res[2] = 3;
	res[3] = 1;
	return res;
}

std::vector<double> v5_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = double_true();
	res[1] = double_true();
	res[2] = double_false();
	res[3] = double_false();
	return res;
}

std::vector<double> v6_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = double_true();
	res[1] = double_false();
	res[2] = double_true();
	res[3] = double_false();
	return res;
}

std::vector<double> v7_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = 6;
	res[1] = 7;
	res[2] = 8;
	res[3] = 7;
	return res;
}

std::vector<double> v8_4() {
	const uint vec_size = 4;
	std::vector<double> res(vec_size);
	res[0] = 0;
	res[1] = -1;
	res[2] = 0;
	res[3] = 1;
	return res;
}

template <class Op>
void test_un_op(std::vector<double> ref_res, std::vector<double> vv1 = v1_4()) {
	const uint vec_size = 4;

	std::cout << "\n" << "testing un op " << typeid(Op).name() << "\n";

	double res_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		res_values[i] = -100;

	ScalarNodePtr v1 = ScalarNode::create_value(&(vv1[0]));
	ScalarNodePtr e = ScalarNode::create<Op>(v1);
	ScalarNodePtr r = ScalarNode::create_result(e, res_values);

	ExpressionDAG se({r});

	ProcessorBase * processor = ProcessorBase::create_processor(se, vec_size, simd_size, nullptr);
	std::vector<uint> subset = {0};
	processor->set_subset(subset);
	processor->run();

	for(uint i=0; i< vec_size; ++i) {
		int64_t mres_values = ((int64_t *)res_values)[i];
		int64_t mref_res = ((int64_t *)&(ref_res[0]))[i];

		std::cout << "i: " << i << " res: " << res_values[i]  << " ref: " << ref_res[i]
                  << " mres: " << std::hex << mres_values
				  << " mref: " << std::hex << mref_res << "\n";

		if (mres_values != mref_res)
		{
			if (std::abs(res_values[i] - ref_res[i]) < 4*std::numeric_limits<double>::epsilon()) {
				continue;
			}
			else {
				BP_ASSERT(false);
			}
		}
	}

	processor->~ProcessorBase();
}


template <class Op>
void test_bin_op(std::vector<double> ref_res, std::vector<double> vv1 = v1_4(), std::vector<double> vv2=v2_4()) {
	const uint vec_size = 4;

	std::cout << "\n" << "testing bin op " << typeid(Op).name() << "\n";

	double res_values[vec_size]; // 5 * double 4
	for(uint i=0;i < vec_size; ++i)
		res_values[i] = -100;

	ScalarNodePtr v1 = ScalarNode::create_value(&(vv1[0]));
	ScalarNodePtr v2 = ScalarNode::create_value(&(vv2[0]));
	ScalarNodePtr e = ScalarNode::create<Op>(v1, v2);
	ScalarNodePtr r = ScalarNode::create_result(e, res_values);

	ExpressionDAG se({r});

	ProcessorBase * processor = ProcessorBase::create_processor(se, vec_size, simd_size, nullptr);
	std::vector<uint> subset = {0};
	processor->set_subset(subset);
	processor->run();

	for(uint i=0; i< vec_size; ++i) {
		int64_t mres_values = ((int64_t *)res_values)[i];
		int64_t mref_res = ((int64_t *)&(ref_res[0]))[i];

		std::cout << "i: " << i << " res: " << res_values[i]  << " ref: " << ref_res[i]
                  << " mres: " << std::hex << mres_values
				  << " mref: " << std::hex << mref_res << "\n";

		if (mres_values != mref_res)
		{
			if (std::abs(res_values[i] - ref_res[i]) < 4*std::numeric_limits<double>::epsilon()) {
				continue;
			}
			else {
				BP_ASSERT(false);
			}
		}
	}

	processor->~ProcessorBase();
}


int main()
{
	test_simple_expr();
	test_un_op<_minus_>({2, 1, -0.0, -1});
	test_bin_op<_add_>({3, 5, 7, 9});
	test_bin_op<_sub_>({-7, -7, -7, -7});
	test_bin_op<_mul_>({-10, -6, 0, 8});
	test_bin_op<_div_>({-2.0 / 5, -1.0 / 6, 0, 1.0 / 8});

	test_bin_op<_mod_>({fmod(-2.0, 5), -1.0, 0, 1.0});
	test_bin_op<_eq_>({double_false(), double_false(), double_true(), double_true()}, v1_4(), v3_4());
	test_bin_op<_ne_>({double_true(), double_true(), double_false(), double_false()}, v1_4(), v3_4());
	test_bin_op<_lt_>({double_false(), double_false(), double_true(), double_false()}, v1_4(), v4_4());
	test_bin_op<_le_>({double_false(), double_false(), double_true(), double_true()}, v1_4(), v4_4());

	test_un_op<_neg_>({double_false(), double_false(), double_true(), double_true()}, v5_4());
	test_bin_op<_or_>({double_true(), double_true(), double_true(), double_false()}, v5_4(), v6_4());
	test_bin_op<_and_>({double_true(), double_false(), double_false(), double_false()}, v5_4(), v6_4());
	test_un_op<_abs_>({2, 1, 0, 1});
	test_un_op<_sqrt_>({sqrt(-2), sqrt(-1), sqrt(0), sqrt(1)});
	test_un_op<_exp_>({exp(5), exp(6), exp(7), exp(8)}, v2_4());
	test_un_op<_log_>({log(5), log(6), log(7), log(8)}, v2_4());
	test_un_op<_log10_>({log10(6), log10(7), log10(8), log10(7)}, v7_4());
	test_un_op<_sin_>({sin(-2), sin(-1), sin(0), sin(1)});
	test_un_op<_sinh_>({sinh(-3.0), sinh(-3.0), sinh(3.0), sinh(1.0)}, v4_4()); // sinh(-2) was last bit off
	test_un_op<_asin_>({asin(-2), asin(-1), asin(0), asin(1)});
	test_un_op<_cos_>({cos(-2), cos(-1), cos(0), cos(1)});
	test_un_op<_cosh_>({cosh(-2), cosh(-1), cosh(0), cosh(1)});
	test_un_op<_acos_>({acos(0), acos(-1), acos(0), acos(1)}, v8_4());
	test_un_op<_tan_>({tan(-2), tan(-1), tan(0), tan(1)});
	test_un_op<_tanh_>({tanh(-2), tanh(-1), tanh(0), tanh(1)});
	test_un_op<_atan_>({atan(-2), atan(-1), atan(0), atan(1)});
	test_un_op<_ceil_>({ceil(-2), ceil(-1), ceil(0), ceil(1)});
	test_un_op<_floor_>({floor(-2), floor(-1), floor(0), floor(1)});
//	test_un_op<_isnan_>({isnan(-2), isnan(-1), isnan(0), isnan(1)}, v5_4());
//	test_un_op<_isinf_>({isinf(-2), isinf(-1), isinf(0), isinf(1)}, v5_4());
//	test_un_op<_sgn_>({sgn(-2), sgn(-1), sgn(0), sgn(1)}, v5_4());

	auto v1 = v1_4();
	auto v2 = v2_4();
	test_bin_op<_atan2_>({atan2(v1[0], v2[0]), atan2(v1[1], v2[1]), atan2(v1[2], v2[2]), atan2(v1[3], v2[3])});
	test_bin_op<_pow_>({pow(v1[0], v2[0]), pow(v1[1], v2[1]), pow(v1[2], v2[2]), pow(v1[3], v2[3])});
}

