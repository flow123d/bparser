# BParser

C++ parsing library for HPC.

Features:
- Python expressions including Numpy like arrays
- compilation to bytecode, arena allocator for maximal memory locality
- expression evaluation on vectors, amortizatiton of the bytecode interpreter 
- usage of SIMD for maximal Flops on modern CPUs 

## Build
SIMD instructions in particular AVX is supported from certain versions of compilers, e.g.
GCC 4.6, Clang 3.5

## Usage:

```
	using namespace bparser;
	
	// Define own value vectors, preferably aligned.
	uint vec_size = 8;
	double v1[vec_size * 3];
	double v2[vec_size * 3];
	double vres[vec_size * 3];

	// Create parser, give the size of the value spaces.
	// That is maximal allocated space. Actual values and 
	// active subset can be changed between evaluations of the expression.
	Parser p(vec_size);
	// parse an expression.
	p.parse("1 * v1 + cs1 * v2");

	// Get symbols used in the expressions.
	std::vector<std::string> symbols = p.symbols();
	std::cout << "Symbols: " << print_vector(p.symbols()) << "\n";
	// Set constants and variables
	
	// "cs1" constant with shape {}, i.e. scalar and values {2}.
	p.set_constant("cs1", {}, 	{2});
	// "cv1" vector constant with shape {3}
	p.set_constant("cv1", {3}, 	{1, 2, 3});
	// "v1" variable with shape {3}; v1 is pointer to the value space
	p.set_variable("v1", {3}, v1);
	p.set_variable("v2", {3}, v2);
	// Set the result variable (the last line of the expression)
	p.set_variable("_result_", {3}, vres);

	// Compile the expression into internall processor.
	p.compile();
	// Set arbitrary subset of the SIMD blocks in the maximal values space.
	// Here the full subset.
	p.set_subset({0, 1});
	// Evaluate
	p.run();
	// Result in the 'vres' value space.
	std::cout << print_vec(vres, 3*vec_size);
```

## Grammar
BParser grammar tries to follow Python 3.6 [expression grammar](https://docs.python.org/3.6/reference/expressions.html). 
In order to perform all array operations statically before compilation  
we restrict indexing and slicing to int literals instead of general expressions.

## Speed

[Intel(R) Core(TM) i7-6820HQ](https://en.wikichip.org/wiki/intel/core_i7/i7-6820hq):

freq. 2.7GHz, turbo boost 3.4GHz, memory bandwith 31.79 GiB/s
L1 data 128KiB, instruction 128KiB,  8 way associative
L2 1MiB, 4 way associative
L3 8MiB
supports AVX instructions, 4 double operations in single tact

theoretical computational power is 10.8 double precision GFLOPs
necessary bandwith for a 2 argument operation (i.e. 12 doubles per SIMD instruction) is
120GiB/s that is about four times what is available so at least 4 operations should be performed on 
the data available in the cache. TODO: experiments with various size of vectors.



## Acknowledgements

- Grammar based on [Henri Menke: boost_matheval](https://github.com/hmenke/boost_matheval)
- Other essential pieces taken from 
[a commented expression parser](https://stackoverflow.com/questions/47354226/how-to-provider-user-with-autocomplete-suggestions-for-given-boostspirit-gramm/47383910#47383910)


## Design

### Vector
Data + subset. 
consists of: SIMD vec, e.g double4 pointer and subset pointer
In debug mode we should store also (maximal) size of the subset and check access to it.
Otherwise the subset size is given by Workspace.

### Workspace
- set of vectors
- buffer for temporary vectors
- buffer for const vectors
- const subset
- single non-const subset
- actual subset size
- all placed in single memory chunk (arena allocator)

### ScalarNode
base class:
- arrity 
- pointers to input ScalarNodes
- index of the node result in the processors workspace (can be shared with some other node)
  !! need to mark somehow if the operation works in place, e.g. +=.
  Then number of inputs is the same as number of eval parameters.
- method add_input to set inputs 
- number of unprocessed dependent nodes (if this is zero we can reuse the temporary storage)

operation specific:
- op code
- store_result = true/false
- call base constructor with appropriate arrity 
- static eval method to perform operation on doubles
  
  
  
  


### Processor - structure for efficient execution of opeartions

- operations are encoded as structs holding an operation code and 
  indices of at most three operands
- table of operands consists of Vectors
- Vector consists of: SIMD vec, e.g double4 pointer and subset pointer

- Processor is created through a static method that takes the ScalarNode tree as the input parameter
  - performs topological sort of the tree (just DFN with a heuristic to reuse the teomporaries slots)
  - count: 
    - input vectors (Vector = data ptr + subset ptr)
    - constants
    - temporaries
    - subset patterns (different pointers)
    - operetions
  - compute and allocate memory
  - construct the Processor, set its pointers to temporaries, constants, subsets
  - create operation list

- Processor is consturucted from a tree of ScalarNodes
- There is specific derived class for every kind of operation.
- common features:
  
### Boost spirit resources
- [Important FAQs](http://boost-spirit.com/home/articles/doc-addendum/faq/#identifiers)
  

  

