

#ifndef INCLUDE_PROCESSOR_BASE_HH_
#define INCLUDE_PROCESSOR_BASE_HH_


#include "arena_alloc.hh"
#include "expression_dag.hh"
#include <vector>

namespace bparser
{

typedef std::shared_ptr<ArenaAlloc> ArenaAllocPtr;
    


struct ProcessorBase {
	virtual void run() = 0;
	virtual void set_subset(std::vector<uint> const &subset) = 0;

	ProcessorBase(ArenaAllocPtr arena)
	: arena_(arena) {
		
	}

	virtual ~ProcessorBase() {
	}

	virtual ArenaAllocPtr get_arena(){
		return arena_;
	}
	
	inline static ProcessorBase *create_processor(details::ExpressionDAG &se, uint vector_size, uint simd_size = 0, ArenaAllocPtr arena = nullptr);

	ArenaAllocPtr arena_;
};



} // namespace bparser

#endif /*INCLUDE_PROCESSOR_BASE_HH_*/