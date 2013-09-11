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

#include "x86/ArrayMathSSE.h"

#if defined(AM_USE_X86) && defined(AM_HAS_SSE)

#include <xmmintrin.h>

#include <cmath>

namespace arraymath {

//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = OP::op(x, *y++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(y) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  __m128 _x = _mm_set_ps1(x);
  if (aligned) {
    for (; length >= 16; length -= 16) {
      _mm_store_ps(dst, OP::opSSE(_x, _mm_load_ps(y)));
      _mm_store_ps(dst + 4, OP::opSSE(_x, _mm_load_ps(y + 4)));
      _mm_store_ps(dst + 8, OP::opSSE(_x, _mm_load_ps(y + 8)));
      _mm_store_ps(dst + 12, OP::opSSE(_x, _mm_load_ps(y + 12)));
      dst += 16; y += 16;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE(_x, _mm_load_ps(y)));
      dst += 4; y += 4;
    }
  }
  else {
    for (; length >= 16; length -= 16) {
      _mm_store_ps(dst, OP::opSSE(_x, _mm_loadu_ps(y)));
      _mm_store_ps(dst + 4, OP::opSSE(_x, _mm_loadu_ps(y + 4)));
      _mm_store_ps(dst + 8, OP::opSSE(_x, _mm_loadu_ps(y + 8)));
      _mm_store_ps(dst + 12, OP::opSSE(_x, _mm_loadu_ps(y + 12)));
      dst += 16; y += 16;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE(_x, _mm_loadu_ps(y)));
      dst += 4; y += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(x, *y++);
  }
}

template <class OP>
void op_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = OP::op(*x++, *y++);
  }

  // Check alignment.
  int alignment = ((reinterpret_cast<size_t>(x) & 15) == 0 ? 1 : 0) |
                  ((reinterpret_cast<size_t>(y) & 15) == 0 ? 2 : 0);

  // 2) Main SSE loop (handle different alignment cases).
  switch (alignment) {
    case 1:
      // x aligned
      for (; length >= 4; length -= 4) {
        _mm_store_ps(dst, OP::opSSE(_mm_load_ps(x), _mm_loadu_ps(y)));
        dst += 4; x += 4; y += 4;
      }
      break;
    case 2:
      // y aligned
      for (; length >= 4; length -= 4) {
        _mm_store_ps(dst, OP::opSSE(_mm_loadu_ps(x), _mm_load_ps(y)));
        dst += 4; x += 4; y += 4;
      }
      break;
    case 3:
      // x aligned && y aligned
      for (; length >= 4; length -= 4) {
        _mm_store_ps(dst, OP::opSSE(_mm_load_ps(x), _mm_load_ps(y)));
        dst += 4; x += 4; y += 4;
      }
      break;
    default:
      // None aligned
      for (; length >= 4; length -= 4) {
        _mm_store_ps(dst, OP::opSSE(_mm_loadu_ps(x), _mm_loadu_ps(y)));
        dst += 4; x += 4; y += 4;
      }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(*x++, *y++);
  }
}


//-----------------------------------------------------------------------------
// Operation implementations.
//-----------------------------------------------------------------------------

struct AddOP {
  static float32 op(float32 a, float32 b) {
    return a + b;
  }
  static __m128 opSSE(__m128 a, __m128 b) {
    return _mm_add_ps(a, b);
  }
};

struct SubOP {
  static float32 op(float32 a, float32 b) {
    return a - b;
  }
  static __m128 opSSE(__m128 a, __m128 b) {
    return _mm_sub_ps(a, b);
  }
};

struct MulOP {
  static float32 op(float32 a, float32 b) {
    return a * b;
  }
  static __m128 opSSE(__m128 a, __m128 b) {
    return _mm_mul_ps(a, b);
  }
};

void ArrayMathSSE::add_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<AddOP>(dst, x, y, length);
}

void ArrayMathSSE::add_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<AddOP>(dst, x, y, length);
}

void ArrayMathSSE::sub_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<SubOP>(dst, x, y, length);
}

void ArrayMathSSE::sub_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<SubOP>(dst, x, y, length);
}

void ArrayMathSSE::mul_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<MulOP>(dst, x, y, length);
}

void ArrayMathSSE::mul_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<MulOP>(dst, x, y, length);
}

void ArrayMathSSE::mulCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main SSE loop.
  __m128 _xr = _mm_set1_ps(xReal);
  __m128 _xi = _mm_set1_ps(xImag);
  for (; length >= 4; length -= 4) {
    __m128 _yr = _mm_loadu_ps(yReal);
    __m128 _yi = _mm_loadu_ps(yImag);
    _mm_storeu_ps(dstReal, _mm_sub_ps(_mm_mul_ps(_xr, _yr), _mm_mul_ps(_xi, _yi)));
    _mm_storeu_ps(dstImag, _mm_add_ps(_mm_mul_ps(_xr, _yi), _mm_mul_ps(_xi, _yr)));
    dstReal += 4; dstImag += 4; yReal += 4; yImag += 4;
  }

  // 2) Tail loop.
  float32 xr = xReal, xi = xImag;
  while (length--) {
    float32 yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathSSE::mulCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main SSE loop.
  for (; length >= 4; length -= 4) {
    __m128 _xr = _mm_loadu_ps(xReal);
    __m128 _xi = _mm_loadu_ps(xImag);
    __m128 _yr = _mm_loadu_ps(yReal);
    __m128 _yi = _mm_loadu_ps(yImag);
    _mm_storeu_ps(dstReal, _mm_sub_ps(_mm_mul_ps(_xr, _yr), _mm_mul_ps(_xi, _yi)));
    _mm_storeu_ps(dstImag, _mm_add_ps(_mm_mul_ps(_xr, _yi), _mm_mul_ps(_xi, _yr)));
    dstReal += 4; dstImag += 4; xReal += 4; xImag += 4; yReal += 4; yImag += 4;
  }

  // 2) Tail loop.
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathSSE::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = std::sqrt(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  const __m128 kZero = _mm_setzero_ps();
  const __m128 kMinus3 = _mm_set_ps1(-3.0f);
  const __m128 kMinus05 = _mm_set_ps1(-0.5f);
  if (aligned) {
    for (; length >= 4; length -= 4) {
      __m128 a = _mm_load_ps(x);
      __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
      a = _mm_mul_ps(r, a);
      r = _mm_mul_ps(a, r);
      r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
      _mm_store_ps(dst, r);
      dst += 4; x += 4;
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      __m128 a = _mm_loadu_ps(x);
      __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
      a = _mm_mul_ps(r, a);
      r = _mm_mul_ps(a, r);
      r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
      _mm_store_ps(dst, r);
      dst += 4; x += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::sqrt(*x++);
  }
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE
