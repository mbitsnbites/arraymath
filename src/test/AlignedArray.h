// -*- Mode: c++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//------------------------------------------------------------------------------
// ArrayMath - an array math library
//------------------------------------------------------------------------------
// Copyright(c) 2013 Marcus Geelnard
//
// This software is provided 'as-is', without any express or implied warranty.
// In no event will the authors be held liable for any damages arising from the
// use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not claim
//    that you wrote the original software. If you use this software in a
//    product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution.
//------------------------------------------------------------------------------

#ifndef _ARRAYMATH_TEST_ALIGNEDARRAY_H
#define _ARRAYMATH_TEST_ALIGNEDARRAY_H

#include <cstddef>

namespace test {

class AlignedArray {
  public:
    AlignedArray(size_t size, size_t alignment) : m_base(0), m_array(0), m_size(0) {
      allocate(size, alignment);
    }

    AlignedArray(size_t size) : m_base(0), m_array(0), m_size(0) {
      allocate(size, kDefaultAlignment);
    }

    ~AlignedArray() {
      deallocate();
    }

    void allocate(size_t size, size_t alignment) {
      deallocate();

      // Allocate the requested array size + additional room for alignment.
      m_base = new float[size + alignment];
      if (!m_base)
        return;
      m_size = size;

      // Start by aligning the buffer to kMaxAlignment.
      size_t adjust = (reinterpret_cast<size_t>(m_base) / sizeof(float)) & (alignment - 1);
      if (adjust)
        adjust = alignment - adjust;
      m_array = m_base + adjust;
    }

    float* get() const {
      return m_array;
    }

    size_t size() const {
      return m_size;
    }

    float& operator[](int i) const {
      return m_array[i];
    }

  private:
    void deallocate() {
      if (m_base) {
        delete[] m_base;
        m_base = 0;
        m_array = 0;
        m_size = 0;
      }
    }

    float* m_base;
    float* m_array;
    size_t m_size;

    AlignedArray& operator=(const AlignedArray&);
    AlignedArray(const AlignedArray&);

    static const size_t kDefaultAlignment = 16;
};

} // namespace test

#endif // _ARRAYMATH_TEST_ALIGNEDARRAY_H
