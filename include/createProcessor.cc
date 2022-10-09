#include "createProcessor.hh"
#include "processorBase.hh"
#include "processor.hh"

namespace bparser{
using namespace details;

inline static uint get_simd_size()
{
	if (__builtin_cpu_supports("avx512f"))
	{
		return 8;
	}
	if (__builtin_cpu_supports("avx2"))
	{
		return 4;
	}
	if (__builtin_cpu_supports("avx"))
	{
		return 4;
	}
	if (__builtin_cpu_supports("sse"))
	{
		return 2;
	}
	else
	{
		return 1;
	}
}

inline ProcessorBase * ProcessorBase::create_processor(details::ExpressionDAG &se, uint vector_size,  uint simd_size, ArenaAllocPtr arena) {
	if (simd_size == 0) {
		simd_size = get_simd_size();
	}
	//std::cout << "n_nodes: " << sorted_nodes.size() << " n_vec: " << se.n_vectors() << "\n";

	switch (simd_size) {
		case 2:
		{
			return create_processor_<Vec2d>(se, vector_size, simd_size, arena);
		} break;
		case 4:
		{
			return create_processor_AVX2(se, vector_size, simd_size, arena);
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

}