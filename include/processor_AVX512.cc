#include "processor.hh"
#include "assert.hh"
#include "arena_alloc.hh"
#include "expression_dag.hh"

#include "VCL_v2_include_AVX512.hh"

namespace bparser{


    template<>
    ProcessorBase * create_processor_<Vec8d>(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena) {
        uint simd_bytes = sizeof(double) * simd_size;
        ExpressionDAG::NodeVec & sorted_nodes = se.sort_nodes();
        uint simd_bytes1 = sizeof(Vec8d);
        // std::cout << simd_bytes1 << "!=" << simd_bytes << "\n";
        BP_ASSERT(simd_bytes1 == simd_bytes);
        uint vec_size = (vector_size / simd_size);
        uint est = 
                align_size(simd_bytes, sizeof(Processor<Vec<Vec8d>>)) +
                align_size(simd_bytes, sizeof(uint) * vector_size) +
                align_size(simd_bytes, se.temp_end * sizeof(Vec<Vec8d>)) +    // temporaries
                align_size(simd_bytes, sizeof(Vec8d) * vec_size * (se.temp_end - se.values_copy_end)) +  // vec_copy, same as temporaries
                align_size(simd_bytes, sizeof(Vec8d) * vec_size * (se.values_copy_end - se.values_end)) + // vector values (probably not neccessary to allocate)
                align_size(simd_bytes, sizeof(Vec8d) * se.constants_end ) +
                align_size(simd_bytes, sizeof(Operation) * (sorted_nodes.size() + 64) );

        // std::cout << "Estimated memory in processor: " << est << std::endl;

        if (arena == nullptr)
            arena = std::make_shared<ArenaAlloc>(simd_bytes, est);
        else
            BP_ASSERT(arena->size_ >= est);
        return arena->create<Processor<Vec<Vec8d>>>(arena, se, vec_size);
    }

        
    ProcessorBase * create_processor_AVX512(ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena) {
        return create_processor_<Vec8d>(se, vector_size, simd_size, arena);
    }


}
#include "createProcessor.hh"