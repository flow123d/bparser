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


size_t align_size(size_t al, size_t size) {
	return (size + al -1) / al * al;
}

struct ArenaAlloc {
	ArenaAlloc(std::size_t alignment, std::size_t size)
	: alignment_(alignment),
	  size_(0)
	{
		size_ = align_size(alignment_, size);
		base_ = (char *)memalign(alignment_, size_);
		ptr_ = base_;
		//std::cout << "arena: " << (void *)base_ << " size: " << size_ << "\n";
	}

	void destroy() {
		free(base_);
	}

	void * allocate(std::size_t size) {
		size = align_size(alignment_, size);
		ASSERT((char *)ptr_ + size <= (char *)base_ + size_);
		void * ptr = ptr_;
		ptr_ += size;
		//std::cout << "allocated: " << ptr << " size: " << size << "\n";
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
