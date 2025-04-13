#include "processor.hh"
#include "expression_dag.hh"
#include "instrset_detect.hh"

namespace bparser{

    inline uint get_simd_size() {
        int i_set = b_instrset_detect();

        if (i_set >= 9)      // min AVX512F
        {
            return 8;
        }
        else if (i_set >= 8) // min AVX2
        {
            return 4;
        }
        else if (i_set >= 5) // min SSE4.1
        {
            return 2;
        }
        else                // no vectorization
        {
            return 1;
        }
    }

    ProcessorBase * ProcessorBase::create_processor(ExpressionDAG &se, uint vector_size, uint simd_size, PatchArenaPtr arena) {
        if (simd_size == 0) {
            simd_size = get_simd_size();
        }

        switch (simd_size) {
            case 2:
            {
                return create_processor_<Vec2d>(se, vector_size, simd_size, arena);
            } break;
            case 4:
            { 
                return create_processor_<Vec4d>(se, vector_size, simd_size, arena);
            } break;
            case 8:
            {
                return create_processor_<Vec8d>(se, vector_size, simd_size, arena);
            } break;
            default:
            {
                return create_processor_<double>(se, vector_size, 1, arena);
            } break;
        }
    }

} // bparser namespace
