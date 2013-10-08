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

//------------------------------------------------------------------------------
// TODO: We get a significant performance hit for the run in/out loops. We
// should implement pure AVX solutions wherever possible!
//------------------------------------------------------------------------------

#include "x86/ArrayMathAVX.h"

#if defined(AM_USE_X86) && defined(AM_HAS_AVX)

#include <immintrin.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace arraymath {

//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = OP::op(x, *y++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(y) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  __m256 _x = _mm256_set1_ps(x);
  if (aligned) {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, OP::opAVX(_x, _mm256_load_ps(y)));
      _mm256_store_ps(dst + 8, OP::opAVX(_x, _mm256_load_ps(y + 8)));
      dst += 16; y += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_x, _mm256_load_ps(y)));
      dst += 8; y += 8;
    }
  }
  else {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, OP::opAVX(_x, _mm256_loadu_ps(y)));
      _mm256_store_ps(dst + 8, OP::opAVX(_x, _mm256_loadu_ps(y + 8)));
      dst += 16; y += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_x, _mm256_loadu_ps(y)));
      dst += 8; y += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(x, *y++);
  }
}

template <class OP>
void op_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = OP::op(*x++, *y++);
  }

  // Check alignment.
  int alignment = ((reinterpret_cast<size_t>(x) & 31) == 0 ? 1 : 0) |
                  ((reinterpret_cast<size_t>(y) & 31) == 0 ? 2 : 0);

  // 2) Main AVX loop (handle different alignment cases).
  switch (alignment) {
    case 1:
      // x aligned
      for (; length >= 8; length -= 8) {
        _mm256_store_ps(dst, OP::opAVX(_mm256_load_ps(x), _mm256_loadu_ps(y)));
        dst += 8; x += 8; y += 8;
      }
      break;
    case 2:
      // y aligned
      for (; length >= 8; length -= 8) {
        _mm256_store_ps(dst, OP::opAVX(_mm256_loadu_ps(x), _mm256_load_ps(y)));
        dst += 8; x += 8; y += 8;
      }
      break;
    case 3:
      // x aligned && y aligned
      for (; length >= 8; length -= 8) {
        _mm256_store_ps(dst, OP::opAVX(_mm256_load_ps(x), _mm256_load_ps(y)));
        dst += 8; x += 8; y += 8;
      }
      break;
    default:
      // None aligned
      for (; length >= 8; length -= 8) {
        _mm256_store_ps(dst, OP::opAVX(_mm256_loadu_ps(x), _mm256_loadu_ps(y)));
        dst += 8; x += 8; y += 8;
      }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(*x++, *y++);
  }
}

template <class OP>
void op_f32_a(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = OP::op(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  if (aligned) {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, OP::opAVX(_mm256_load_ps(x)));
      _mm256_store_ps(dst + 8, OP::opAVX(_mm256_load_ps(x + 8)));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_mm256_load_ps(x)));
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, OP::opAVX(_mm256_loadu_ps(x)));
      _mm256_store_ps(dst + 8, OP::opAVX(_mm256_loadu_ps(x + 8)));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_mm256_loadu_ps(x)));
      dst += 8; x += 8;
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
  static __m256 opAVX(__m256 a, __m256 b) {
    return _mm256_add_ps(a, b);
  }
};

struct SubOP {
  static float32 op(float32 a, float32 b) {
    return a - b;
  }
  static __m256 opAVX(__m256 a, __m256 b) {
    return _mm256_sub_ps(a, b);
  }
};

struct MulOP {
  static float32 op(float32 a, float32 b) {
    return a * b;
  }
  static __m256 opAVX(__m256 a, __m256 b) {
    return _mm256_mul_ps(a, b);
  }
};

struct DivOP {
  static float32 op(float32 a, float32 b) {
    return a / b;
  }
  static __m256 opAVX(__m256 a, __m256 b) {
    return _mm256_div_ps(a, b);
  }
};

struct FloorOP {
  static float32 op(float32 a) {
    return std::floor(a);
  }
  static __m256 opAVX(__m256 a) {
    return _mm256_floor_ps(a);
  }
};

struct CeilOP {
  static float32 op(float32 a) {
    return std::ceil(a);
  }
  static __m256 opAVX(__m256 a) {
    return _mm256_ceil_ps(a);
  }
};

struct RoundOP {
  static float32 op(float32 a) {
    return std::floor(a + 0.5f);
  }
  static __m256 opAVX(__m256 a) {
    return _mm256_round_ps(a, _MM_FROUND_NINT);
  }
};

struct FractOP {
  static float32 op(float32 a) {
    return a - std::floor(a);
  }
  static __m256 opAVX(__m256 a) {
    return _mm256_sub_ps(a, _mm256_floor_ps(a));
  }
};

void ArrayMathAVX::add_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<AddOP>(dst, x, y, length);
}

void ArrayMathAVX::add_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<AddOP>(dst, x, y, length);
}

void ArrayMathAVX::sub_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<SubOP>(dst, x, y, length);
}

void ArrayMathAVX::sub_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<SubOP>(dst, x, y, length);
}

void ArrayMathAVX::mul_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<MulOP>(dst, x, y, length);
}

void ArrayMathAVX::mul_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<MulOP>(dst, x, y, length);
}

void ArrayMathAVX::mulCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main AVX loop.
  __m256 _xr = _mm256_set1_ps(xReal);
  __m256 _xi = _mm256_set1_ps(xImag);
  for (; length >= 8; length -= 8) {
    __m256 _yr = _mm256_loadu_ps(yReal);
    __m256 _yi = _mm256_loadu_ps(yImag);
    _mm256_storeu_ps(dstReal, _mm256_sub_ps(_mm256_mul_ps(_xr, _yr), _mm256_mul_ps(_xi, _yi)));
    _mm256_storeu_ps(dstImag, _mm256_add_ps(_mm256_mul_ps(_xr, _yi), _mm256_mul_ps(_xi, _yr)));
    dstReal += 8; dstImag += 8; yReal += 8; yImag += 8;
  }

  // 2) Tail loop.
  float32 xr = xReal, xi = xImag;
  while (length--) {
    float32 yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathAVX::mulCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main AVX loop.
  for (; length >= 8; length -= 8) {
    __m256 _xr = _mm256_loadu_ps(xReal);
    __m256 _xi = _mm256_loadu_ps(xImag);
    __m256 _yr = _mm256_loadu_ps(yReal);
    __m256 _yi = _mm256_loadu_ps(yImag);
    _mm256_storeu_ps(dstReal, _mm256_sub_ps(_mm256_mul_ps(_xr, _yr), _mm256_mul_ps(_xi, _yi)));
    _mm256_storeu_ps(dstImag, _mm256_add_ps(_mm256_mul_ps(_xr, _yi), _mm256_mul_ps(_xi, _yr)));
    dstReal += 8; dstImag += 8; xReal += 8; xImag += 8; yReal += 8; yImag += 8;
  }

  // 2) Tail loop.
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathAVX::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<DivOP>(dst, x, y, length);
}

void ArrayMathAVX::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<DivOP>(dst, x, y, length);
}

void ArrayMathAVX::madd_f32_saa(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = x * *y++ + *z++;
  }

  // Check alignment.
  int alignment = ((reinterpret_cast<size_t>(y) & 31) == 0 ? 1 : 0) |
                  ((reinterpret_cast<size_t>(z) & 31) == 0 ? 2 : 0);

  // 2) Main AVX loop (handle different alignment cases).
  __m256 _x = _mm256_set1_ps(x);
  switch (alignment) {
    case 1:
      // y aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_x, _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; y += 8; z += 8;
      }
      break;
    case 2:
      // z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_x, _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; y += 8; z += 8;
      }
      break;
    case 3:
      // y aligned && z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_x, _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; y += 8; z += 8;
      }
      break;
    default:
      // None aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_x, _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; y += 8; z += 8;
      }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = x * *y++ + *z++;
  }
}

void ArrayMathAVX::madd_f32_aaa(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = *x++ * *y++ + *z++;
  }

  // Check alignment.
  // Check alignment.
  int alignment = ((reinterpret_cast<size_t>(x) & 31) == 0 ? 1 : 0) |
                  ((reinterpret_cast<size_t>(y) & 31) == 0 ? 2 : 0) |
                  ((reinterpret_cast<size_t>(z) & 31) == 0 ? 4 : 0);

  // 2) Main AVX loop (handle different alignment cases).
  switch (alignment) {
    case 1:
      // x aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 2:
      // y aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 3:
      // x aligned && y aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 4:
      // z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 5:
      // x aligned && z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 6:
      // y aligned && z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    case 7:
      // x aligned && y aligned && z aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_load_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
      break;
    default:
      // None aligned
      for (; length >= 8; length -= 8) {
        __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_loadu_ps(y));
        _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
        dst += 8; x += 8; y += 8; z += 8;
      }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = *x++ * *y++ + *z++;
  }
}

void ArrayMathAVX::abs_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = std::abs(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main AVX loop.
  static const __m256 kMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
  if (aligned) {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, _mm256_and_ps(_mm256_load_ps(x), kMask));
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, _mm256_and_ps(_mm256_loadu_ps(x), kMask));
      dst += 8; x += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::abs(*x++);
  }
}

void ArrayMathAVX::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = std::sqrt(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  static const __m256 kZero = _mm256_setzero_ps();
  static const __m256 kMinus3 = _mm256_set1_ps(-3.0f);
  static const __m256 kMinus05 = _mm256_set1_ps(-0.5f);
  if (aligned) {
    for (; length >= 8; length -= 8) {
      __m256 a = _mm256_load_ps(x);
      __m256 r = _mm256_and_ps(_mm256_rsqrt_ps(a), _mm256_cmp_ps(kZero, a, _CMP_NEQ_UQ));
      a = _mm256_mul_ps(r, a);
      r = _mm256_mul_ps(a, r);
      r = _mm256_mul_ps(_mm256_mul_ps(kMinus05, a), _mm256_add_ps(kMinus3, r));
      _mm256_store_ps(dst, r);
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      __m256 a = _mm256_loadu_ps(x);
      __m256 r = _mm256_and_ps(_mm256_rsqrt_ps(a), _mm256_cmp_ps(kZero, a, _CMP_NEQ_UQ));
      a = _mm256_mul_ps(r, a);
      r = _mm256_mul_ps(a, r);
      r = _mm256_mul_ps(_mm256_mul_ps(kMinus05, a), _mm256_add_ps(kMinus3, r));
      _mm256_store_ps(dst, r);
      dst += 8; x += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = std::sqrt(*x++);
  }
}

void ArrayMathAVX::ceil_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CeilOP>(dst, x, length);
}

void ArrayMathAVX::floor_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FloorOP>(dst, x, length);
}

float32 ArrayMathAVX::max_f32(const float32 *x, size_t length) {
  float32 result = -std::numeric_limits<float>::infinity();

  // 1) Align x to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(x) & 31) && length--) {
    result = std::max(result, *x++);
  }

  // 2) Main AVX loop.
  __m256 _result = _mm256_set1_ps(result);
  __m256 _result2 = _result;
  for (; length >= 16; length -= 16) {
    _result = _mm256_max_ps(_result, _mm256_load_ps(x));
    _result2 = _mm256_max_ps(_result2, _mm256_load_ps(x + 8));
    x += 16;
  }
  _result = _mm256_max_ps(_result, _result2);
  for (; length >= 8; length -= 8) {
    _result = _mm256_max_ps(_result, _mm256_load_ps(x));
    x += 8;
  }

  // 3) Horizontal max of SIMD register.
  union {
    __m256 v;
    float32 f[8];
  } u;
  u.v = _result;
  result = u.f[0];
  for (int k = 1; k < 8; ++k) {
    result = std::max(result, u.f[k]);
  }

  // 4) Tail loop.
  while (length--) {
    result = std::max(result, *x++);
  }

  return result;
}

float32 ArrayMathAVX::min_f32(const float32 *x, size_t length) {
  float32 result = std::numeric_limits<float>::infinity();

  // 1) Align x to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(x) & 31) && length--) {
    result = std::min(result, *x++);
  }

  // 2) Main AVX loop.
  __m256 _result = _mm256_set1_ps(result);
  __m256 _result2 = _result;
  for (; length >= 16; length -= 16) {
    _result = _mm256_min_ps(_result, _mm256_load_ps(x));
    _result2 = _mm256_min_ps(_result2, _mm256_load_ps(x + 8));
    x += 16;
  }
  _result = _mm256_min_ps(_result, _result2);
  for (; length >= 8; length -= 8) {
    _result = _mm256_min_ps(_result, _mm256_load_ps(x));
    x += 8;
  }

  // 3) Horizontal min of SIMD register.
  union {
    __m256 v;
    float32 f[8];
  } u;
  u.v = _result;
  result = u.f[0];
  for (int k = 1; k < 8; ++k) {
    result = std::min(result, u.f[k]);
  }

  // 4) Tail loop.
  while (length--) {
    result = std::min(result, *x++);
  }

  return result;
}

void ArrayMathAVX::round_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<RoundOP>(dst, x, length);
}

void ArrayMathAVX::clamp_f32(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  __m256 _xMin = _mm256_set1_ps(xMin);
  __m256 _xMax = _mm256_set1_ps(xMax);
  if (aligned) {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_load_ps(x))));
      _mm256_store_ps(dst + 8, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_load_ps(x + 8))));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_load_ps(x))));
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 16; length -= 16) {
      _mm256_store_ps(dst, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_loadu_ps(x))));
      _mm256_store_ps(dst + 8, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_loadu_ps(x + 8))));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, _mm256_max_ps(_xMin, _mm256_min_ps(_xMax, _mm256_loadu_ps(x))));
      dst += 8; x += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }
}

void ArrayMathAVX::fract_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FractOP>(dst, x, length);
}

void ArrayMathAVX::fill_f32(float32 *dst, float32 value, size_t length) {
  if (length >= 16) {
    // 1) Align dst to a 32-byte boundary.
    size_t num_unaligned = (reinterpret_cast<size_t>(dst) & 31) >> 2;
    num_unaligned = num_unaligned ? 8 - num_unaligned : 0;
    length -= num_unaligned;
    while (num_unaligned--) {
      *dst++ = value;
    }

    // 2) Main AVX loop.
    __m256 _value = _mm256_set1_ps(value);
    for (; length >= 32; length -= 32) {
      _mm256_store_ps(dst, _value);
      _mm256_store_ps(dst + 8, _value);
      _mm256_store_ps(dst + 16, _value);
      _mm256_store_ps(dst + 24, _value);
      dst += 32;
    }
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, _value);
      dst += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = value;
  }
}

void ArrayMathAVX::ramp_f32(float32 *dst, float32 first, float32 last, size_t length) {
  if (length == 0) {
    return;
  }
  if (length == 1) {
    *dst = first;
    return;
  }

  float32 step = (last - first) / static_cast<float32>(length - 1);
  float32 k = 0.0f;

  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }

  // 2) Main AVX loop.
  static const __m256 kEight = _mm256_set1_ps(8.0f);
  static const __m256 kIdxRamp = _mm256_set_ps(7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f);
  __m256 _first = _mm256_set1_ps(first);
  __m256 _step = _mm256_set1_ps(step);
  __m256 _k = _mm256_add_ps(_mm256_set1_ps(k), kIdxRamp);
  size_t mainLoopSize = (length / 8) * 8;
  for (; length >= 8; length -= 8) {
    _mm256_store_ps(dst, _mm256_add_ps(_first, _mm256_mul_ps(_step, _k)));
    _k = _mm256_add_ps(_k, kEight);
    dst += 8;
  }
  k += mainLoopSize;

  // 3) Tail loop.
  while (length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }
}

void ArrayMathAVX::sign_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = *x++ < 0.0f ? -1.0f : 1.0f;
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  static const __m256 kZero = _mm256_setzero_ps();
  static const __m256 kOne = _mm256_set1_ps(1.0f);
  static const __m256 kMinusOne = _mm256_set1_ps(-1.0f);
  if (aligned) {
    for (; length >= 16; length -= 16) {
      __m256 mask = _mm256_cmp_ps(kZero, _mm256_load_ps(x), _CMP_LT_OS);
      _mm256_store_ps(dst, _mm256_blendv_ps(kMinusOne, kOne, mask));
      mask = _mm256_cmp_ps(kZero, _mm256_load_ps(x + 8), _CMP_LT_OS);
      _mm256_store_ps(dst + 8, _mm256_blendv_ps(kMinusOne, kOne, mask));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      __m256 mask = _mm256_cmp_ps(kZero, _mm256_load_ps(x), _CMP_LT_OS);
      _mm256_store_ps(dst, _mm256_blendv_ps(kMinusOne, kOne, mask));
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 16; length -= 16) {
      __m256 mask = _mm256_cmp_ps(kZero, _mm256_loadu_ps(x), _CMP_LT_OS);
      _mm256_store_ps(dst, _mm256_blendv_ps(kMinusOne, kOne, mask));
      mask = _mm256_cmp_ps(kZero, _mm256_loadu_ps(x + 8), _CMP_LT_OS);
      _mm256_store_ps(dst + 8, _mm256_blendv_ps(kMinusOne, kOne, mask));
      dst += 16; x += 16;
    }
    for (; length >= 8; length -= 8) {
      __m256 mask = _mm256_cmp_ps(kZero, _mm256_loadu_ps(x), _CMP_LT_OS);
      _mm256_store_ps(dst, _mm256_blendv_ps(kMinusOne, kOne, mask));
      dst += 8; x += 8;
    }
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = *x++ < 0.0f ? -1.0f : 1.0f;
  }
}

void ArrayMathAVX::sampleLinear_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  __m256 _maxIdx = _mm256_set1_ps(static_cast<float32>(maxIdx));
  const __m256 kZero = _mm256_set1_ps(0.0f);
  const __m256 kOne = _mm256_set1_ps(1.0f);

  union {
    __m256i v;
    int i[8];
  } idx1, idx2;

  union {
    __m256 v;
    float f[8];
  } p1, p2;

  // 1) Main AVX loop.
  for (; length >= 8; length -= 8) {
    __m256 _t2 = _mm256_loadu_ps(t);
    _t2 = _mm256_max_ps(kZero, _mm256_min_ps(_maxIdx, _t2));
    __m256 _w = _mm256_floor_ps(_t2);

    idx1.v = _mm256_cvtps_epi32(_w);
    idx2.v = _mm256_cvtps_epi32(_mm256_min_ps(_maxIdx, _mm256_add_ps(_w, kOne)));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 8; ++k) {
      p1.f[k] = x[idx1.i[k]];
      p2.f[k] = x[idx2.i[k]];
    }

    _w = _mm256_sub_ps(_t2, _w);
    _mm256_storeu_ps(dst, _mm256_add_ps(p1.v, _mm256_mul_ps(_w, _mm256_sub_ps(p2.v, p1.v))));

    dst += 8; t += 8;
  }

  // 2) Tail loop.
  while (length--) {
    float32 t2 = *t++;
    t2 = t2 < 0 ? 0 : t2 > maxIdx ? maxIdx : t2;
    size_t idx = std::floor(t2);
    float32 w = t2 - static_cast<float32>(idx);
    float32 p1 = x[idx];
    float32 p2 = x[idx < maxIdx ? idx + 1 : maxIdx];
    *dst++ = p1 + w * (p2 - p1);
  }
}

void ArrayMathAVX::sampleLinearRepeat_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32 xLengthF = static_cast<float32>(xLength);
  float32 xLengthFInv = 1.0f / xLengthF;
  __m256 _maxIdx = _mm256_set1_ps(static_cast<float32>(maxIdx));
  __m256 _xLengthF = _mm256_set1_ps(xLengthF);
  __m256 _xLengthFInv = _mm256_set1_ps(xLengthFInv);
  const __m256 kZero = _mm256_set1_ps(0.0f);
  const __m256 kOne = _mm256_set1_ps(1.0f);

  union {
    __m256i v;
    int i[8];
  } idx1, idx2;

  union {
    __m256 v;
    float f[8];
  } p1, p2;

  // 1) Main AVX loop.
  for (; length >= 8; length -= 8) {
    __m256 _t2 = _mm256_loadu_ps(t);
    _t2 = _mm256_sub_ps(_t2, _mm256_mul_ps(_mm256_floor_ps(_mm256_mul_ps(_t2, _xLengthFInv)), _xLengthF));
    __m256 _w = _mm256_floor_ps(_t2);

    idx1.v = _mm256_cvtps_epi32(_w);
    idx2.v = _mm256_cvtps_epi32(_mm256_and_ps(_mm256_cmp_ps(_w, _maxIdx, _CMP_NEQ_UQ), _mm256_add_ps(_w, kOne)));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 8; ++k) {
      p1.f[k] = x[idx1.i[k]];
      p2.f[k] = x[idx2.i[k]];
    }

    _w = _mm256_sub_ps(_t2, _w);
    _mm256_storeu_ps(dst, _mm256_add_ps(p1.v, _mm256_mul_ps(_w, _mm256_sub_ps(p2.v, p1.v))));

    dst += 8; t += 8;
  }

  // 2) Tail loop.
  while (length--) {
    float32 t2 = *t++;
    t2 = t2 - std::floor(t2 * xLengthFInv) * xLengthF;
    size_t idx = std::floor(t2);
    float32 w = t2 - static_cast<float32>(idx);
    float32 p1 = x[idx];
    float32 p2 = x[idx < maxIdx ? idx + 1 : 0];
    *dst++ = p1 + w * (p2 - p1);
  }
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_AVX
