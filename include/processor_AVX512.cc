#include "processor.hh"
#include "VCL_v2_include_AVX512.hh"

namespace bparser{

    template<>
    ProcessorBase * create_processor_<Vec8d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
}
