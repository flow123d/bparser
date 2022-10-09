#include "processorBase.hh"
#include <stdint.h>
#include "expression_dag.hh"

namespace bparser{

    inline static uint get_simd_size();

    inline ProcessorBase * create_processor(details::ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena);


}