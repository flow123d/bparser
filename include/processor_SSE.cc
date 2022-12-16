#include "processor.hh"
#include "VCL_v2_include_SSE.hh"

namespace bparser{

    template<>
    ProcessorBase * create_processor_<Vec2d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
 }
 