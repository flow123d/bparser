#include "processor.hh"

// save diagnostic state
#pragma GCC diagnostic push 

// turn off the specific warning.
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"

#undef MAX_VECTOR_SIZE
#define MAX_VECTOR_SIZE 256

#include "vectorclass.h"

// turn the warnings back on
#pragma GCC diagnostic pop

namespace bparser{

    template<>
    ProcessorBase * create_processor_<Vec4d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
}
