#include <stdio.h>
#include <chrono>
#include <iostream>
#include "..//third_party/VCL_v2/vectorclass.h"
#include "..//third_party/VCL_v2/instrset_detect.cpp"

template <typename VecType>
void eval(VecType &res, VecType a, VecType b) 
{
	res = a < b;
}
//#define MAX_VECTOR_SIZE 256
int main() 
{
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < 1000000; i++)
    {
        // Vec8d dsd(1.2, 2.4, 3.6, 4.8, 5.6, 4.2, 23.4, 234.6);
        // Vec8d res8dsd = dsd + dsd;

        // Vec4d dsd1(1.2, 2.4, 3.6, 4.8);
        // Vec4d dsd2(5.6, 4.2, 23.4, 234.6);
        // Vec4d res4dsd1 = dsd1 + dsd1;
        // Vec4d res4dsd2 = dsd2 + dsd2;
    }
    
    Vec4d aa(0.0, 1.0, 20.0, 30.0);
    Vec4d bb(10.0, 20.0, 3.0, 4.0);
    Vec4d rr;

    Vec4d cc;
    Vec4d &cc_ref = cc;
    cc_ref = aa < bb;

    eval<Vec4d>(rr, aa, bb);

    printf("\nResult rr = ");
    for(int i = 0; i < rr.size(); i++)
    {
        printf("%f ", rr[i]);
    }
    printf("\nResult cc_ref = ");
    for(int i = 0; i < cc_ref.size(); i++)
    {
        printf("%f ", cc_ref[i]);
    }

    Vec4db dd = aa < bb;
    cc_ref = dd;

    printf("\nResult cc_ref = ");
    for(int i = 0; i < cc_ref.size(); i++)
    {
        printf("%f ", cc_ref[i]);
    }

    return 0;

    auto end_time = std::chrono::high_resolution_clock::now();

	double time  = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    printf("%f \n", time);

    return 0;

    typedef std::conditional<INSTRSET >= 7, Vec4d, Vec2d>::type tmp;
    typedef std::conditional<INSTRSET >= 9, Vec8d, tmp>::type VecType;
    
    VecType res;

    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < 1000000; i++)
    {

    //VecType test(1.2, 2.4, 3.6, 4.8);
    //VecType trst(5.6, 4.2, 23.4, 234.6);

    //res = test + trst;
    }
    end_time = std::chrono::high_resolution_clock::now();

	time  = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    printf("%f \n", time);


    int is = instrset_detect();
    int instrs = INSTRSET;
    printf("%i ", is);
    printf("%i ", instrs);
    printf("%i ", INSTRSET);
    printf("%i ", MAX_VECTOR_SIZE);
    printf("%li ", sizeof(Vec4d));
    printf("%li ", sizeof(Vec4db));
    printf("\nResult res_test = ");
    for(int i = 0; i < res.size(); i++)
    {
        printf("%f ", res[i]);
    }

    return 0;
    // vec of 16 integers
    Vec16i a(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    Vec16i b(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);

    // vecs of 4
    Vec4i c(1, 2, 3, 4);
    Vec4d d(1.2, 2.4, 3.6, 4.8);
    Vec4f f(1.2, 2.4, 3.6, 4.8);

    Vec16i res16i = a * b;
    Vec4i res4i = 5 + c;
    Vec4d res4d = pow(d, 3);
    Vec4f res4f = sqrt(f);

    printf("\nResult res16i = ");
    for(int i = 0; i < res16i.size(); i++)
    {
        printf("%i ", res16i[i]);
    }

    printf("\nResult res4i = ");
    for(int i = 0; i < res4i.size(); i++)
    {
        printf("%i ", res4i[i]);
    }

    printf("\nResult res4d = ");
    for(int i = 0; i < res4d.size(); i++)
    {
        printf("%f ", res4d[i]);
    }

    printf("\nResult res4f = ");
    for(int i = 0; i < res4f.size(); i++)
    {
        printf("%f ", res4f[i]);
    }

    return 0;
}
//*/

/*
template<typename VEC_T> //udělat na Vec4d

void printVector(VEC_T & v, const char * prefix) 
{
    std::cout << prefix << " ";
    for(int i = 0; i < VEC_T::length(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

int main()
{
    // vec's of 4 doubles
    Vec4d a(12, 13, 14, 15);
    Vec4d b(1, 2, 3, 4);

    Vec4d c;
    Vec4d d;
    Vec4d e;

   // Vec4d c(2, 4, 8, 16);
   // Vec4d d(1.2, 2.4, 3.6, 4.8);
   // Vec4d e(10.2, 20.4, 30.6, 40.8);
    

// Labels are useful for searching through assembly listings.
test1_start:
    c = a + b;
test1_end:
    printVector(c, "c:" );

test2_start:
    d = a + 2 * b;
test2_end:
    printVector(d, "d:" );

test3_start:
    e = a * b;
test3_end:
    printVector(e, "e:" );

    return 0;
}
*/


///*