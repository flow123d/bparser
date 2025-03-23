/*
 * aligned_alloc.hh
 *
 *  Created on: Mar 23, 2025
 *      Author: LV
 */

/*
* Abstracts the memory allocation and deallocation functions for crosscompile reasons (Windows / Linux)
*/

#ifndef INCLUDE_ALIGNED_ALLOC_HH
#define INCLUDE_ALIGNED_ALLOC_HH


#include <cstdlib>
#include <malloc.h>

namespace bparser {
	//https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170&viewFallbackFrom=vs-2019#note_M
	inline void* align_alloc(size_t size, size_t alignment) {
#if defined(_WIN32)
		return _aligned_malloc(size, alignment);
#else
		return memalign(alignment, size);
#endif
	}

	inline void align_free(void* p) {
#if defined(_WIN32)
		_aligned_free(p);
#else
		free(p);
#endif
	}
}


#endif // !INCLUDE_ALIGNED_ALLOC_HH
