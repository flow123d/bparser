/*
 * linear_allocator.hh
 *
 *  Created on: Dec 31, 2019
 *      Author: jb
 */

#ifndef INCLUDE_LINEAR_ALLOCATOR_HH_
#define INCLUDE_LINEAR_ALLOCATOR_HH_

#include <cstddef>
#include <cassert>
#include <new>
#include "aligned_mallocations.hpp"


template <typename T>
class LinearAllocator
{
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    //using propagate_on_container_copy_assignment = std::true_type;
    //using propagate_on_container_move_assignment = std::true_type;
    //using propagate_on_container_swap = std::true_type;

    LinearAllocator(std::size_t count = 64)
        : m_memUsed(0),
        m_memStartAddress(nullptr)
    {
        allocate(count);
    }
    ~LinearAllocator()
    {
        clear();
    }

    template <class U>
    LinearAllocator(const LinearAllocator<U>&) noexcept
    {}

    /// \brief allocates memory equal to # count objects of type T
    pointer allocate(std::size_t count)
    {
        if (count > std::size_t(-1) / sizeof(T))
        {
            throw std::bad_alloc{};
        }
        if (m_memStartAddress != nullptr)
        {
            alignedFree(m_memStartAddress);
        }
        m_memUsed = count * sizeof(T);
        m_memStartAddress = static_cast<pointer>(alignedMalloc(m_memUsed, alignof(T)));
        return m_memStartAddress;
    }
    /// \brief deallocates previously allocated memory
    /// \brief Linear/arena allocators do not support free() operations. Use clear() instead.
    void deallocate([[maybe_unused]] pointer p, [[maybe_unused]] std::size_t count) noexcept
    {
        //assert(false);
        clear();
    }

    /// \brief simply resets memory
    void clear()
    {
        if (m_memStartAddress != nullptr)
        {
            alignedFree(m_memStartAddress);
            m_memStartAddress = nullptr;
        }
        this->m_memUsed = 0;
    }

    /// \brief GETTERS
    pointer getStartAddress() const
    {
        return this->m_memStartAddress;
    }
    std::size_t getUsedMemory() const
    {
        return this->m_memUsed;
    }
private:
    std::size_t m_memUsed;
    pointer m_memStartAddress;
};

template <class T, class U>
bool operator==(const LinearAllocator<T> &, const LinearAllocator<U> &)
{
    return true;
}

template <class T, class U>
bool operator!=(const LinearAllocator<T> &, const LinearAllocator<U> &)
{
    return false;
}




#endif /* INCLUDE_LINEAR_ALLOCATOR_HH_ */
