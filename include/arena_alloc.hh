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
#include "aligned_alloc.hh"
#include "arena_resource.hh"

namespace bparser {


inline size_t align_size(size_t al, size_t size) {
	return (size + al -1) / al * al;
}

struct ArenaAlloc {
	

	//Creates a wrapper of PatchArena for backwards compatibility with BParser
	ArenaAlloc(PatchArena& existing_arena) : arena(&existing_arena),buffer(nullptr) {
		;
	}
	//Creates a wrapper with a new PatchArena with the specified memory alignment and size
	//However AssemblyArena might be the correct class to create 
	ArenaAlloc(std::size_t alignment, std::size_t size)
	//: alignment_(alignment),
	//  size_(0)
	: size_(align_size(alignment, size)) //We cannot access this->arena->buffer_size_. However it *should* not change and is currently only used by one assert. Maybe create getter? -LV
	{
		buffer = align_alloc(alignment, size_);
		arena = new PatchArena(buffer, size_, alignment);
		/*size_ = align_size(alignment_, size);
		base_ = (char*)align_alloc(alignment_, size_);
		BP_ASSERT(base_ != nullptr);
		ptr_ = base_;
		//std::cout << "arena begin: " << (void *)base_ << " end: " << end() << std::endl;
		*/
	}
	
	~ArenaAlloc() {
        destroy();
    }

	inline void destroy() {
		//align_free(base_);
		if (buffer != nullptr) {
			align_free(buffer);
		}
	}

	/*void* end() {
		return base_ + size_;
	}*/

	inline void* allocate(std::size_t size) {
		/*
		size = align_size(alignment_, size);
		void * ptr = ptr_;
		ptr_ += size;
		BP_ASSERT(ptr_ <= end());
		//std::cout << "allocated: " << ptr << " end: " << (void *)ptr_ << " aend: " << end() << "\n";
		return ptr;
		*/
		return arena->allocate(size);
		
	}

	template <class T, typename... Args>
	T* create(Args&&... args) {
		
		void * ptr = allocate(sizeof(T));
		return new (ptr) T(std::forward<Args>(args)...);
		
	}

	template <class T>
	T* create_array(uint n_items) {
		/*
		void * ptr = allocate(sizeof(T) * n_items);
		return new (ptr) T[n_items];
		*/
		return arena->allocate_simd<T>(n_items);
	}

	inline std::size_t get_size() {
		return size_; //arena->buffer_size would be more appropriate
	}
	
	//std::size_t alignment_;
	//std::size_t size_;
	//char * base_;
	//char * ptr_;
protected:
	PatchArena* arena;
	void* buffer;
	std::size_t size_;
};

} // namespace bparser


#endif /* INCLUDE_ARENA_ALLOC_HH_ */
