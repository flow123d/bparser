// vec size and p names must stay the same!

#ifndef NITPICK_IDE_IGNORE
# include "nitpick_include.hh"
#endif //This block will not be compiled during generation/run. But it helps the IDE understand the code

constexpr uint vec_size = 8;
constexpr uint max_vec_size = vec_size;

//MEM_ALLOC(name, 3*vec_size, 2) //2,2,2,2,2,2...
//MEM_ALLOC_INDEX(namei, 3 * vec_size) // 0,1,2,3,4,5...
//MEM_ALLOC_LINEAR(namel, 3 * vec_size) // 1,2,3,4,5,6...
constexpr uint vres_size = 3 * vec_size;
MEM_ALLOC(vres, vres_size, NAN)

Parser p(max_vec_size);

//P_SET_CONSTANT(name, ARG({2,2}), ARG({1,2,3,4})) //2x2 matrix
//P_SET_VARIABLE(v1, ARG({3}), namei) //vector from namei double*
P_SET_VARIABLE(_result_, { }, vres);

p.parse("1+1");