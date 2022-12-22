#include "processor.hh"

//test if gnuc
#ifdef __GNUC__

// save diagnostic state
#pragma GCC diagnostic push 

// turn off the specific warning.
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

// end if
#endif

#undef MAX_VECTOR_SIZE
#define MAX_VECTOR_SIZE 512

#include "vectorclass.h"

//test if gnuc
#ifdef __GNUC__

// turn the warnings back on
#pragma GCC diagnostic pop

// end if
#endif


namespace bparser{

    template<>
    ProcessorBase * create_processor_<Vec8d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
}
