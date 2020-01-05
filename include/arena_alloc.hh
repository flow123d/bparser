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
namespace details {

size_t align_size(size_t al, size_t size) {
	return (size + al -1) / al * al;
}

struct ArenaAlloc {
	ArenaAlloc(std::size_t align_size, std::size_t size)
	: alignment_(align_size),
	  size_(size)
	{
		base_ = memalign(alignment_, size_);
		ptr_ = base_;
	}

	void destroy() {
		free(base_);
	}

	void * allocate(std::size_t size) {
		size = align_size(alignment_, alignment_);
		ASSERT(ptr_ + size <= base_ + size_);
		void * ptr = ptr_;
		ptr_ += size;
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
	void * base_;
	void * ptr_;
};

} // namespace details
} // namespace bparser



#endif /* INCLUDE_ARENA_ALLOC_HH_ */
