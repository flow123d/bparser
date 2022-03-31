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

## Expression syntax
BParser grammar tries to follow Python 3.6 [expression grammar](https://docs.python.org/3.6/reference/expressions.html). 

Whole program consists of sequence of assignements separated by ';' finished with the result expression.

```
    a = var + 1;
    b = a + var; c=2*a + b;
    a + b + c
```

White space doesn't metter. Input variables (as the `var`) are detected and must be provided through the 
`set_variable` or `set_constant` methods (see 'Usage' section).
New variables are defined by the assignment. Identifiers are formed from a-z,A-Z,_,0-9 
charactes and must not start with a digit.

Supported unary operators: 
    - `not` : boolean inversion
    - `+` : unary plus
    - `-` : unary minus
    
Supported operators (sorted from the highest precedence):
    - `**` : power
    - unary `+`, unary `-`
    - `*`, `/`, `//`, `%`, `@` : multiply, division, integer division, modulo, matrix multiplication
    - `+`, `-`
    - `<`, `>`, `<=`, `>=`, `==`, `!=`
    - `not`
    - `and`
    - `or`
    - `<true_expr> if <condition> else <false_expr>` : conditional expression
    
Chained comparison is suported, e. g.

`a < b < c == d`

is evaluated as:

`a<b and b<c and c==d`

Operands can be: numbers, variables, function calls, and array subscriptions.

Supported functions:
   - `abs`, `ceil`,  `floor`, `sgn`
   - `rad2deg` (radians to degrees), `deg2rad` (degrees to radians),
   - `acos`, `asin`, `atan`, `atan2`, `cos`, `sin`, `tan`,
   - `isinf` (equal to `+/-inf`, `isnan` (equal to `+/- Nan),
   - `log`, `log10`, `log2`, `exp`, `cosh`, `sinh`
   - `sqrt`, `power`
   - `eye`, `zeros`, `ones`, `full` : array construction
   - `flatten`
   - `minimum`, `maximum`
   
Arrays: vectors, matrices, tensors.

Construction:
`[1,2,3]` : vector, shape {3}
`[[1,2], [3,4]]` : matrix, shape {2,3}
 
Subscription, slices:
`a[0, 1]` : element on row 0, column 1 of a matrix
`a[[5, 0, 2]]` : subvector of elements at indices 5, 0, 2
`a[-1]` : negative indices are treated as `axis_size +  index`

`a[0:10:2]` : eleemnts at indices 0,2,4,6,8
`a[-1:-2:-1]` : last two elements in reversed order

`a[None, 1]` : vector of shape `{n}` broadcasted to a matrix of shape `{1,n}`


## Syntax grammar
program:
    (assignment)+ expression 

assignment:
    identifier '=' expression ';'

expression:
    | disjunction 'if' disjunction 'else' expression
    | disjunction

disjunction:
    | conjunction ('or' conjunction)+
    | conjunction

conjunction:
    | inversion ('and' inversion)+
    | inversion
    
inversion:
    | 'not' inversion 
    | comparison

comparison:
    | bitwise_or compare_op_bitwise_or_pair+ 
    | bitwise_or

compare_op_bitwise_or_pair:
    cmp_op sum

cmp_op:
    {'==', '!=', '<', '<=', '>', '>='}

sum:
    | sum '+' term 
    | sum '-' term 
    | term
    
term:
    | term '*' factor 
    | term '/' factor 
    | term '//' factor 
    | term '%' factor 
    | term '@' factor 
    | factor
    
factor:
    | '+' factor 
    | '-' factor 
    | power
    
power:
    | primary '**' factor 
    | primary

primary:
    | primary '.' NAME 
    | primary genexp 
    | primary '(' [arguments] ')' 
    | primary '[' slices ']' 
    | atom

slices:
    | slice !',' 
    | ','.slice+ [',']
    
slice:
    | [expression] ':' [expression] [':' [expression] ] 
    | named_expression 
    
atom:
    | NAME
    | 'True' 
    | 'False' 
    | 'None' 
    | NUMBER

strings: STRING+ 
list:
    | '[' [star_named_expressions] ']'
Literals: 
  - integer: e.g. 1
  - float: exponential notation, e.g. 1.0; 1.0e-0
  - bool: 'True' or 'False'
  - variable name: any string formed from a-z,A-Z,_,0-9, not starting with a digit
  

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

- Boost spirit grammar based on [Henri Menke: boost_matheval](https://github.com/hmenke/boost_matheval)
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
  

  

