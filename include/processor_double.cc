#include "processor.hh"

namespace bparser{

    template<>
    ProcessorBase * create_processor_<double>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);
}
