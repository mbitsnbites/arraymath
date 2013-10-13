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

#include "x86/ArrayMathSSE4.h"

#if defined(AM_USE_X86) && defined(AM_HAS_SSE4)

#include <smmintrin.h>

#include <cmath>

namespace arraymath {

namespace {

//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_a(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = OP::op(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 15) == 0;

  // 2) Main SSE4 loop (handle different alignment cases).
  if (aligned) {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE4(_mm_load_ps(x)));
      dst += 4; x += 4;
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE4(_mm_loadu_ps(x)));
      dst += 4; x += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(*x++);
  }
}


//-----------------------------------------------------------------------------
// Operation implementations.
//-----------------------------------------------------------------------------

struct FloorOP {
  static float32 op(float32 a) {
    return std::floor(a);
  }
  static __m128 opSSE4(__m128 a) {
    return _mm_round_ps(a, _MM_FROUND_FLOOR);
  }
};

struct CeilOP {
  static float32 op(float32 a) {
    return std::ceil(a);
  }
  static __m128 opSSE4(__m128 a) {
    return _mm_round_ps(a, _MM_FROUND_CEIL);
  }
};

struct RoundOP {
  static float32 op(float32 a) {
    return std::floor(a + 0.5f);
  }
  static __m128 opSSE4(__m128 a) {
    return _mm_round_ps(a, _MM_FROUND_NINT);
  }
};

} // anonymous namespace


//-----------------------------------------------------------------------------
// Class methods.
//-----------------------------------------------------------------------------

void ArrayMathSSE4::floor_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FloorOP>(dst, x, length);
}

void ArrayMathSSE4::ceil_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CeilOP>(dst, x, length);
}

void ArrayMathSSE4::round_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<RoundOP>(dst, x, length);
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE4
