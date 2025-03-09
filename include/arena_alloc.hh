/*
 * arena_alloc.hh
 *
 *  Created on: Jan 5, 2020
 *      Author: jb
 */

#ifndef INCLUDE_ARENA_ALLOC_HH_
#define INCLUDE_ARENA_ALLOC_HH_

#include <cstdlib>
#include <utility>
#include <malloc.h>

namespace bparser {


inline size_t align_size(size_t al, size_t size) {
	return (size + al -1) / al * al;
}

struct ArenaAlloc {
	ArenaAlloc(std::size_t alignment, std::size_t size)
	: alignment_(alignment),
	  size_(0)
	{
		size_ = align_size(alignment_, size);
		base_ = (char *)_aligned_malloc(size_, alignment_); //https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170&viewFallbackFrom=vs-2019#note_M
		BP_ASSERT(base_ != nullptr);
		ptr_ = base_;
		//std::cout << "arena begin: " << (void *)base_ << " end: " << end() << std::endl;
	}
	
	~ArenaAlloc() {
        destroy();
    }

	void destroy() {
		_aligned_free(base_);
	}

	void *end() {
		return base_ + size_;
	}

	void * allocate(std::size_t size) {
		size = align_size(alignment_, size);
		void * ptr = ptr_;
		ptr_ += size;
		BP_ASSERT(ptr_ <= end());
		//std::cout << "allocated: " << ptr << " end: " << (void *)ptr_ << " aend: " << end() << "\n";
		return ptr;
	}

	template <class T, typename... Args>
	T * create(Args&&... args) {
		void * ptr = allocate(sizeof(T));
		return new (ptr) T(std::forward<Args>(args)...);
	}

	template <class T>
	T * create_array(uint n_items) {
		void * ptr = allocate(sizeof(T) * n_items);
		return new (ptr) T[n_items];
	}


	std::size_t alignment_;
	std::size_t size_;
	char * base_;
	char * ptr_;
};

} // namespace bparser


#endif /* INCLUDE_ARENA_ALLOC_HH_ */
