

## Thanks
Grammar based on [Henri Menke: boost_matheval](https://github.com/hmenke/boost_matheval)

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
  

  

  

