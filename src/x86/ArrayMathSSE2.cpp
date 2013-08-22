// -*- Mode: c++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//------------------------------------------------------------------------------
// ArrayMath - an array math library
//------------------------------------------------------------------------------
// Copyright (c) 2013 Marcus Geelnard
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

#include "ArrayMathSSE2.h"

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)

#include <emmintrin.h>

namespace arraymath {

//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  // 1) Align y to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(y) & 15) && length--) {
    *dst++ = OP::op(x, *y++);
  }
  const __m128 *_y = reinterpret_cast<const __m128*>(y);

  // Check alignment.
  bool dstAligned = (reinterpret_cast<size_t>(dst) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  __m128 _x = _mm_set_ps1(x);
  if (dstAligned) {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE2(_x, *_y++));
      dst += 4;
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      _mm_storeu_ps(dst, OP::opSSE2(_x, *_y++));
      dst += 4;
    }
  }

  // 3) Tail loop.
  y = reinterpret_cast<const float32*>(_y);
  while (length--) {
    *dst++ = OP::op(x, *y++);
  }
}

template <class OP>
void op_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  // 1) Align x to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(x) & 15) && length--) {
    *dst++ = OP::op(*x++, *y++);
  }
  const __m128 *_x = reinterpret_cast<const __m128*>(x);

  // Check alignment.
  bool dstAligned = (reinterpret_cast<size_t>(dst) & 15) == 0;
  bool yAligned = (reinterpret_cast<size_t>(y) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  if (dstAligned && yAligned) {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE2(*_x++, _mm_load_ps(y)));
      dst += 4;
      y += 4;
    }
  }
  else if (dstAligned) {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE2(*_x++, _mm_loadu_ps(y)));
      dst += 4;
      y += 4;
    }
  }
  else if (yAligned) {
    for (; length >= 4; length -= 4) {
      _mm_storeu_ps(dst, OP::opSSE2(*_x++, _mm_load_ps(y)));
      dst += 4;
      y += 4;
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      _mm_storeu_ps(dst, OP::opSSE2(*_x++, _mm_loadu_ps(y)));
      dst += 4;
      y += 4;
    }
  }

  // 3) Tail loop.
  x = reinterpret_cast<const float32*>(_x);
  while (length--) {
    *dst++ = OP::op(*x++, *y++);
  }
}


//-----------------------------------------------------------------------------
// Operation implementations.
//-----------------------------------------------------------------------------

struct DivOP {
  static float32 op(float32 a, float32 b) {
    return a / b;
  }
  static __m128 opSSE2(__m128 a, __m128 b) {
    return _mm_div_ps(a, b);
  }
};

void ArrayMathSSE2::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<DivOP>(dst, x, y, length);
}

void ArrayMathSSE2::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<DivOP>(dst, x, y, length);
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE2
