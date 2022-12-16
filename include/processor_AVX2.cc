#include "processor.hh"
#include "VCL_v2_include_AVX2.hh"

namespace bparser{

    template<>
    ProcessorBase * create_processor_<Vec4d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
}
