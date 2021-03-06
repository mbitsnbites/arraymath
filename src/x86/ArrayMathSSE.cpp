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

#include <algorithm>
#include <cmath>
#include <limits>

namespace arraymath {

namespace {

//-----------------------------------------------------------------------------
// Helper functions.
//-----------------------------------------------------------------------------

#include "x86/sse_mathfun.h"

// Reciprocal with 23-bits of precision (_mm_rcp_ps only gives 12 bits of
// precision).
__m128 rcp_23bit(__m128 x) {
  __m128 r = _mm_rcp_ps(x);
  __m128 mask = _mm_cmpneq_ps(x, _mm_setzero_ps());
  __m128 filtered = _mm_and_ps(mask, _mm_mul_ps(_mm_mul_ps(r, r), x));
  return _mm_sub_ps(_mm_add_ps(r, r), filtered);
}


//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
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
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
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

template <class OP>
void op_f32_a(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = OP::op(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  if (aligned) {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE(_mm_load_ps(x)));
      dst += 4; x += 4;
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, OP::opSSE(_mm_loadu_ps(x)));
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

struct DivOP {
  static float32 op(float32 a, float32 b) {
    return a / b;
  }
  static __m128 opSSE(__m128 a, __m128 b) {
    return _mm_mul_ps(a, rcp_23bit(b));
  }
};

struct SinOP {
  static float32 op(float32 a) {
    return std::sin(a);
  }
  static __m128 opSSE(__m128 a) {
    return sin_ps(a);
  }
};

struct CosOP {
  static float32 op(float32 a) {
    return std::cos(a);
  }
  static __m128 opSSE(__m128 a) {
    return cos_ps(a);
  }
};

struct TanOP {
  static float32 op(float32 a) {
    return std::tan(a);
  }
  static __m128 opSSE(__m128 a) {
    __m128 s, c;
    sincos_ps(a, &s, &c);
    return _mm_mul_ps(s, rcp_23bit(c));
  }
};

struct ExpOP {
  static float32 op(float32 a) {
    return std::exp(a);
  }
  static __m128 opSSE(__m128 a) {
    return exp_ps(a);
  }
};

struct LogOP {
  static float32 op(float32 a) {
    return std::log(a);
  }
  static __m128 opSSE(__m128 a) {
    return log_ps(a);
  }
};

} // anonymous namespace


//-----------------------------------------------------------------------------
// Class methods.
//-----------------------------------------------------------------------------

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

void ArrayMathSSE::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<DivOP>(dst, x, y, length);
}

void ArrayMathSSE::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<DivOP>(dst, x, y, length);
}

void ArrayMathSSE::divCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main AVX loop.
  __m128 _xr = _mm_set1_ps(xReal);
  __m128 _xi = _mm_set1_ps(xImag);
  for (; length >= 4; length -= 4) {
    __m128 _yr = _mm_loadu_ps(yReal);
    __m128 _yi = _mm_loadu_ps(yImag);
    __m128 _denom = rcp_23bit(_mm_add_ps(_mm_mul_ps(_yr, _yr), _mm_mul_ps(_yi, _yi)));
    _mm_storeu_ps(dstReal, _mm_mul_ps(_mm_add_ps(_mm_mul_ps(_xr, _yr), _mm_mul_ps(_xi, _yi)), _denom));
    _mm_storeu_ps(dstImag, _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_xi, _yr), _mm_mul_ps(_xr, _yi)), _denom));
    dstReal += 4; dstImag += 4; yReal += 4; yImag += 4;
  }

  // 2) Tail loop.
  float32 xr = xReal, xi = xImag;
  while (length--) {
    float32 yr = *yReal++, yi = *yImag++;
    float32 denom = 1.0f / (yr * yr + yi * yi);
    *dstReal++ = (xr * yr + xi * yi) * denom;
    *dstImag++ = (xi * yr - xr * yi) * denom;
  }
}

void ArrayMathSSE::divCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main AVX loop.
  for (; length >= 4; length -= 4) {
    __m128 _xr = _mm_loadu_ps(xReal);
    __m128 _xi = _mm_loadu_ps(xImag);
    __m128 _yr = _mm_loadu_ps(yReal);
    __m128 _yi = _mm_loadu_ps(yImag);
    __m128 _denom = rcp_23bit(_mm_add_ps(_mm_mul_ps(_yr, _yr), _mm_mul_ps(_yi, _yi)));
    _mm_storeu_ps(dstReal, _mm_mul_ps(_mm_add_ps(_mm_mul_ps(_xr, _yr), _mm_mul_ps(_xi, _yi)), _denom));
    _mm_storeu_ps(dstImag, _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_xi, _yr), _mm_mul_ps(_xr, _yi)), _denom));
    dstReal += 4; dstImag += 4; xReal += 4; xImag += 4; yReal += 4; yImag += 4;
  }

  // 2) Tail loop.
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    float32 denom = 1.0f / (yr * yr + yi * yi);
    *dstReal++ = (xr * yr + xi * yi) * denom;
    *dstImag++ = (xi * yr - xr * yi) * denom;
  }
}

float32 ArrayMathSSE::max_f32(const float32 *x, size_t length) {
  float32 result = -std::numeric_limits<float>::infinity();

  // 1) Align x to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(x) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    result = std::max(result, *x++);
  }

  // 2) Main SSE loop.
  __m128 _result = _mm_set1_ps(result);
  __m128 _result2 = _result;
  __m128 _result3 = _result;
  __m128 _result4 = _result;
  for (; length >= 16; length -= 16) {
    _result = _mm_max_ps(_result, _mm_load_ps(x));
    _result2 = _mm_max_ps(_result2, _mm_load_ps(x + 4));
    _result3 = _mm_max_ps(_result3, _mm_load_ps(x + 8));
    _result4 = _mm_max_ps(_result4, _mm_load_ps(x + 12));
    x += 16;
  }
  _result3 = _mm_max_ps(_result3, _result4);
  _result2 = _mm_max_ps(_result2, _result3);
  _result = _mm_max_ps(_result, _result2);
  for (; length >= 4; length -= 4) {
    _result = _mm_max_ps(_result, _mm_load_ps(x));
    x += 4;
  }

  // 3) Horizontal max of SIMD register.
  union {
    __m128 v;
    float32 f[4];
  } u;
  u.v = _result;
  result = u.f[0];
  for (int k = 1; k < 4; ++k) {
    result = std::max(result, u.f[k]);
  }

  // 4) Tail loop.
  while (length--) {
    result = std::max(result, *x++);
  }

  return result;
}

float32 ArrayMathSSE::min_f32(const float32 *x, size_t length) {
  float32 result = std::numeric_limits<float>::infinity();

  // 1) Align x to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(x) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    result = std::min(result, *x++);
  }

  // 2) Main SSE loop.
  __m128 _result = _mm_set1_ps(result);
  __m128 _result2 = _result;
  __m128 _result3 = _result;
  __m128 _result4 = _result;
  for (; length >= 16; length -= 16) {
    _result = _mm_min_ps(_result, _mm_load_ps(x));
    _result2 = _mm_min_ps(_result2, _mm_load_ps(x + 4));
    _result3 = _mm_min_ps(_result3, _mm_load_ps(x + 8));
    _result4 = _mm_min_ps(_result4, _mm_load_ps(x + 12));
    x += 16;
  }
  _result3 = _mm_min_ps(_result3, _result4);
  _result2 = _mm_min_ps(_result2, _result3);
  _result = _mm_min_ps(_result, _result2);
  for (; length >= 4; length -= 4) {
    _result = _mm_min_ps(_result, _mm_load_ps(x));
    x += 4;
  }

  // 3) Horizontal min of SIMD register.
  union {
    __m128 v;
    float32 f[4];
  } u;
  u.v = _result;
  result = u.f[0];
  for (int k = 1; k < 4; ++k) {
    result = std::min(result, u.f[k]);
  }

  // 4) Tail loop.
  while (length--) {
    result = std::min(result, *x++);
  }

  return result;
}

void ArrayMathSSE::abs_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = std::abs(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 15) == 0;

  // 2) Main SSE loop.
  static const __m128 kMask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
  if (aligned) {
    for (; length >= 8; length -= 8) {
      _mm_store_ps(dst, _mm_and_ps(_mm_load_ps(x), kMask));
      _mm_store_ps(dst + 4, _mm_and_ps(_mm_load_ps(x + 4), kMask));
      dst += 8; x += 8;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, _mm_and_ps(_mm_load_ps(x), kMask));
      dst += 4; x += 4;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm_store_ps(dst, _mm_and_ps(_mm_loadu_ps(x), kMask));
      _mm_store_ps(dst + 4, _mm_and_ps(_mm_loadu_ps(x + 4), kMask));
      dst += 8; x += 8;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, _mm_and_ps(_mm_loadu_ps(x), kMask));
      dst += 4; x += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::abs(*x++);
  }
}

void ArrayMathSSE::absCplx_f32(float32 *dst, const float32 *xReal, const float32 *xImag, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    float32 re = *xReal++;
    float32 im = *xImag++;
    *dst++ = std::sqrt(re * re + im * im);
  }

  // Check alignment.
  int alignment = ((reinterpret_cast<size_t>(xReal) & 15) == 0 ? 1 : 0) |
                  ((reinterpret_cast<size_t>(xImag) & 15) == 0 ? 2 : 0);

  // 2) Main SSE loop (handle different alignment cases).
  static const __m128 kZero = _mm_setzero_ps();
  static const __m128 kMinus3 = _mm_set1_ps(-3.0f);
  static const __m128 kMinus05 = _mm_set1_ps(-0.5f);
  switch (alignment) {
    case 1:
      // xReal aligned
      for (; length >= 4; length -= 4) {
        __m128 re = _mm_load_ps(xReal);
        __m128 im = _mm_loadu_ps(xImag);
        __m128 a = _mm_add_ps(_mm_mul_ps(re, re), _mm_mul_ps(im, im));
        __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
        a = _mm_mul_ps(r, a);
        r = _mm_mul_ps(a, r);
        r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
        _mm_store_ps(dst, r);
        dst += 4; xReal += 4; xImag += 4;
      }
      break;
    case 2:
      // xImag aligned
      for (; length >= 4; length -= 4) {
        __m128 re = _mm_loadu_ps(xReal);
        __m128 im = _mm_load_ps(xImag);
        __m128 a = _mm_add_ps(_mm_mul_ps(re, re), _mm_mul_ps(im, im));
        __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
        a = _mm_mul_ps(r, a);
        r = _mm_mul_ps(a, r);
        r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
        _mm_store_ps(dst, r);
        dst += 4; xReal += 4; xImag += 4;
      }
      break;
    case 3:
      // xReal aligned && xImag aligned
      for (; length >= 4; length -= 4) {
        __m128 re = _mm_load_ps(xReal);
        __m128 im = _mm_load_ps(xImag);
        __m128 a = _mm_add_ps(_mm_mul_ps(re, re), _mm_mul_ps(im, im));
        __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
        a = _mm_mul_ps(r, a);
        r = _mm_mul_ps(a, r);
        r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
        _mm_store_ps(dst, r);
        dst += 4; xReal += 4; xImag += 4;
      }
      break;
    default:
      // None aligned
      for (; length >= 4; length -= 4) {
        __m128 re = _mm_loadu_ps(xReal);
        __m128 im = _mm_loadu_ps(xImag);
        __m128 a = _mm_add_ps(_mm_mul_ps(re, re), _mm_mul_ps(im, im));
        __m128 r = _mm_and_ps(_mm_rsqrt_ps(a), _mm_cmpneq_ps(kZero, a));
        a = _mm_mul_ps(r, a);
        r = _mm_mul_ps(a, r);
        r = _mm_mul_ps(_mm_mul_ps(kMinus05, a), _mm_add_ps(kMinus3, r));
        _mm_store_ps(dst, r);
        dst += 4; xReal += 4; xImag += 4;
      }
  }

  // 3) Tail loop.
  while (length--) {
    float32 re = *xReal++;
    float32 im = *xImag++;
    *dst++ = std::sqrt(re * re + im * im);
  }
}

void ArrayMathSSE::sin_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<SinOP>(dst, x, length);
}

void ArrayMathSSE::cos_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CosOP>(dst, x, length);
}

void ArrayMathSSE::tan_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<TanOP>(dst, x, length);
}

void ArrayMathSSE::exp_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<ExpOP>(dst, x, length);
}

void ArrayMathSSE::log_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<LogOP>(dst, x, length);
}

void ArrayMathSSE::pow_f32_as(float32 *dst, const float32 *x, float32 y, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = std::pow(*x++, y);
  }

  // 2) Main SSE loop.
  // Note: x^y = exp(y*log(x))
  // TODO(m): This can probably be optimized even further.
  __m128 _y = _mm_set1_ps(y);
  for (; length >= 4; length -= 4) {
    __m128 _x = _mm_loadu_ps(x);
    _mm_store_ps(dst, exp_ps(_mm_mul_ps(_y, log_ps(_x))));
    dst += 4; x += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::pow(*x++, y);
  }
}

void ArrayMathSSE::pow_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = std::pow(*x++, *y++);
  }

  // 2) Main SSE loop.
  // Note: x^y = exp(y*log(x))
  // TODO(m): This can probably be optimized even further.
  for (; length >= 4; length -= 4) {
    __m128 _x = _mm_loadu_ps(x);
    __m128 _y = _mm_loadu_ps(y);
    _mm_store_ps(dst, exp_ps(_mm_mul_ps(_y, log_ps(_x))));
    dst += 4; x += 4; y += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::pow(*x++, *y++);
  }
}

void ArrayMathSSE::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
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

void ArrayMathSSE::clamp_f32(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 15) == 0;

  // 2) Main SSE loop (handle different alignment cases).
  __m128 _xMin = _mm_set1_ps(xMin);
  __m128 _xMax = _mm_set1_ps(xMax);
  if (aligned) {
    for (; length >= 8; length -= 8) {
      _mm_store_ps(dst, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_load_ps(x))));
      _mm_store_ps(dst + 4, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_load_ps(x + 4))));
      dst += 8; x += 8;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_load_ps(x))));
      dst += 4; x += 4;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm_store_ps(dst, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_loadu_ps(x))));
      _mm_store_ps(dst + 4, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_loadu_ps(x + 4))));
      dst += 8; x += 8;
    }
    for (; length >= 4; length -= 4) {
      _mm_store_ps(dst, _mm_max_ps(_xMin, _mm_min_ps(_xMax, _mm_loadu_ps(x))));
      dst += 4; x += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }
}

void ArrayMathSSE::fill_f32(float32 *dst, float32 value, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = value;
  }

  // 2) Main SSE loop.
  __m128 _value = _mm_set1_ps(value);
  for (; length >= 16; length -= 16) {
    _mm_store_ps(dst, _value);
    _mm_store_ps(dst + 4, _value);
    _mm_store_ps(dst + 8, _value);
    _mm_store_ps(dst + 12, _value);
    dst += 16;
  }
  for (; length >= 4; length -= 4) {
    _mm_store_ps(dst, _value);
    dst += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = value;
  }
}

void ArrayMathSSE::ramp_f32(float32 *dst, float32 first, float32 last, size_t length) {
  if (AM_UNLIKELY(length == 0)) {
    return;
  }
  if (AM_UNLIKELY(length == 1)) {
    *dst = first;
    return;
  }

  float32 step = (last - first) / static_cast<float32>(length - 1);
  float32 k = 0.0f;

  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }

  // 2) Main SSE loop.
  static const __m128 kFour = _mm_set1_ps(4.0f);
  static const __m128 kIdxRamp = _mm_set_ps(3.0f, 2.0f, 1.0f, 0.0f);
  __m128 _first = _mm_set1_ps(first);
  __m128 _step = _mm_set1_ps(step);
  __m128 _k = _mm_add_ps(_mm_set1_ps(k), kIdxRamp);
  size_t mainLoopSize = (length / 4) * 4;
  for (; length >= 4; length -= 4) {
    _mm_store_ps(dst, _mm_add_ps(_first, _mm_mul_ps(_step, _k)));
    _k = _mm_add_ps(_k, kFour);
    dst += 4;
  }
  k += mainLoopSize;

  // 3) Tail loop.
  while (length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE
