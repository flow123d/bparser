#include <stdio.h>
#include <chrono>
#include <iostream>
#include "test_tools.hh"
#include "VCL_v2_include.hh"
#include "..//third_party/VCL_v2/instrset_detect.cpp"

const int64_t true_value = 0xFFFFFFFFFFFFFFFFLL; //-1LL;
const int64_t false_value = 0x0000000000000000LL;

// template<typename VecType>
// void print_VCL_vector(const VecType & v, const char * prefix)
// {
//     bool first = true;
//     std::cout << prefix << "(";
//     for(int i = 0; i < VecType::size(); i++)
//     {
//         if (first)
//         {
//             std::cout << v[i];
//             first = false;
//             continue;
//         }

//         std::cout << " ; " << v[i];
//     }
//     std::cout << ")" << std::endl;
// }

template <typename VecType>
void eval(VecType &res, VecType a, VecType b)
{
	// as_bool(res) = a < b;
    // res = as_double(a < b);
    res = a - b * truncate(a / b);
}

template<typename bool_type> struct b_to_d;

template<>
struct b_to_d<int64_t> {
    typedef double double_type;
};

template<>
struct b_to_d<Vec2db> {
    typedef Vec2d double_type;
};

template<>
struct b_to_d<Vec4db> {
    typedef Vec4d double_type;
};

template<>
struct b_to_d<Vec8db> {
    typedef Vec8d double_type;
};

template<typename bool_type> union b_to_d_mask;

template<>
union b_to_d_mask<int64_t> {
	int64_t	mask;
	double  value;
};

template<>
union b_to_d_mask<Vec2db> {
	Vec2db	mask;
	Vec2d  value;
};

template<>
union b_to_d_mask<Vec4db> {
	Vec4db	mask;
	Vec4d  value;
};

template<>
union b_to_d_mask<Vec8db> {
	Vec8db	mask;
	Vec8d  value;
};

template<typename bool_type>
inline typename b_to_d<bool_type>::double_type as_double(bool_type in) {
    b_to_d_mask<bool_type> x = {in};
    return x.value;
}

template<typename double_type> struct d_to_b;

template<>
struct d_to_b<double> {
    typedef int64_t bool_type;
};

template<>
struct d_to_b<Vec2d> {
    typedef Vec2db bool_type;
};

template<>
struct d_to_b<Vec4d> {
    typedef Vec4db bool_type;
};

template<>
struct d_to_b<Vec8d> {
    typedef Vec8db bool_type;
};

template<typename double_type> union d_to_b_mask;

template<>
union d_to_b_mask<double> {
	double  value;
	int64_t	mask;
};

template<>
union d_to_b_mask<Vec2d> {
	Vec2d  value;
	Vec2db	mask;
};

template<>
union d_to_b_mask<Vec4d> {
	Vec4d  value;
	Vec4db	mask;
};

template<>
union d_to_b_mask<Vec8d> {
	Vec8d  value;
	Vec8db	mask;
};


template<typename double_type>
inline typename d_to_b<double_type>::bool_type as_bool(double_type in) {
    d_to_b_mask<double_type> x = {in};
    return x.mask;
}


//#define MAX_VECTOR_SIZE 256
int main()
{
    Vec2d a(0.0, 5.0);
    Vec2d b(1.0, 4.0);
    Vec2d r;
    Vec2db x;
    
    x = a > b;
    r = as_double(x);

    print_VCL_vector<Vec2db>(x, "x");
    print_VCL_vector<Vec2d>(r, "r");

    r = select(x, a, b);
    print_VCL_vector<Vec2d>(r, "r");

///////////////////////////////////////////////////////////////
    Vec4d aa(0.0, 1.0, 20.0, 30.0);
    Vec4d bb(10.0, 20.0, 3.0, 4.0);
    Vec4d rr;

    Vec4d cc;
    Vec4d &cc_ref = cc;
    cc_ref = aa < bb;

    eval<Vec4d>(rr, aa, bb);

    print_VCL_vector<Vec4d>(rr, "rr");
    print_VCL_vector<Vec4db>(as_bool(rr), "rr");
    print_VCL_vector<Vec4d>(cc_ref, "cc_ref");

    Vec4db xx = aa < bb;
    cc_ref = xx;

    print_VCL_vector<Vec4d>(cc_ref, "cc_ref");

    cc_ref = as_double(xx);

    print_VCL_vector<Vec4d>(cc_ref, "cc_ref");

    rr = select(xx, aa, bb);
    print_VCL_vector<Vec4d>(rr, "rr");


    std::cout << std::endl;
    std::cout << std::endl;
    Vec4d t_a(3.0, 3.0, 3.0, 3.0);
    Vec4d t_b;
    Vec4d t_c(4.0, 4.0, 4.0, 4.0);
    Vec4db t_bin;

    t_bin = t_a < t_c;
    print_VCL_vector<Vec4db>(t_bin, "t_bin");

    t_b = as_double(t_a < t_c);
    // t_b = true_value;

    print_VCL_vector<Vec4d>(t_a, "t_a");
    print_VCL_vector<Vec4d>(t_b, "t_b");
    print_VCL_vector<d_to_b<Vec4d>::bool_type>(as_bool(t_b), "t_b_bool");
    print_VCL_vector<Vec4d>(t_c, "t_c");

    Vec4d t_r = select(as_bool(t_b), t_a, t_c);
    
    print_VCL_vector<Vec4d>(t_r, "select_res");

    std::cout << std::endl;
    std::cout << std::endl;

////////////////////////////////////////////////////////////////////////
    Vec8d aaa(0.0, 1.0, 20.0, 30.0, 0.0, 1.0, 20.0, 30.0);
    Vec8d bbb(10.0, 20.0, 3.0, 4.0, 0.0, 1.0, 20.0, 30.0);
    Vec8d rrr;
    Vec8db xxx;

    xxx = aaa < bbb;

    print_VCL_vector<Vec8db>(xxx, "xxx");

    rrr = as_double(xxx);

    print_VCL_vector<Vec8d>(rrr, "rrr");

    as_bool(rrr) = xxx;

    print_VCL_vector<d_to_b<Vec8d>::bool_type>(as_bool(rrr), "rrr");

    Vec8d ccc;
    Vec8d &ccc_ref = ccc;
    ccc_ref = as_double(aaa < bbb);

    print_VCL_vector<Vec8d>(ccc_ref, "ccc_ref");
    print_VCL_vector<Vec8db>(as_bool(ccc_ref), "ccc_ref");

    eval<Vec8d>(rrr, aaa, bbb);

    print_VCL_vector<Vec8d>(rrr, "rrr");

    rrr = select(xxx, aaa, bbb);
    print_VCL_vector<Vec8d>(rrr, "rrr");

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < 1000000; i++)
    {
        eval<Vec2d>(r, a, b);
        eval<Vec4d>(rr, aa, bb);
        eval<Vec8d>(rrr, aaa, bbb);   
    }

    auto end_time = std::chrono::high_resolution_clock::now();

	double time  = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    printf("time elapsed: %f \n", time);

    return 0;

    // typedef std::conditional<INSTRSET >= 7, Vec4d, Vec2d>::type tmp;
    // typedef std::conditional<INSTRSET >= 9, Vec8d, tmp>::type VecType;

    // VecType res;

    // start_time = std::chrono::high_resolution_clock::now();
    // for (int i = 1; i < 1000000; i++)
    // {
    //     VecType test(1.2, 2.4, 3.6, 4.8);
    //     VecType trst(5.6, 4.2, 23.4, 234.6);

    //     res = test + trst;
    // }
    // end_time = std::chrono::high_resolution_clock::now();

	// time  = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    // printf("%f \n", time);


    // int is = instrset_detect();
    // int instrs = INSTRSET;
    // printf("%i ", is);
    // printf("%i ", instrs);
    // printf("%i ", INSTRSET);
    // printf("%i ", MAX_VECTOR_SIZE);
    // printf("%li ", sizeof(Vec4d));
    // printf("%li ", sizeof(Vec4db));
    // printf("\nResult res_test = ");
    // for(int i = 0; i < res.size(); i++)
    // {
    //     printf("%f ", res[i]);
    // }

    // return 0;
}


// // Labels are useful for searching through assembly listings.
// test1_start:
//     c = a + b;
// test1_end:
//     print_VCL_vector(c, "c:" );

// test2_start:
//     d = a + 2 * b;
// test2_end:
//     print_VCL_vector(d, "d:" );

// test3_start:
//     e = a * b;
// test3_end:
//     print_VCL_vector(e, "e:" );

