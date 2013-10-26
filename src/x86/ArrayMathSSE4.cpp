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

#include "x86/ArrayMathSSE.h"

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

struct FractOP {
  static float32 op(float32 a) {
    return a - std::floor(a);
  }
  static __m128 opSSE4(__m128 a) {
    return _mm_sub_ps(a, _mm_round_ps(a, _MM_FROUND_FLOOR));
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

void ArrayMathSSE4::fract_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FractOP>(dst, x, length);
}

void ArrayMathSSE4::sampleLinear_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (AM_UNLIKELY(xLength == 0)) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    ArrayMathSSE::fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  __m128 _maxIdx = _mm_set1_ps(static_cast<float32>(maxIdx));
  const __m128 kZero = _mm_set1_ps(0.0f);
  const __m128 kOne = _mm_set1_ps(1.0f);

  union {
    __m128i v;
    int i[4];
  } idx1, idx2;

  union {
    __m128 v;
    float f[4];
  } p1, p2;

  // 1) Main SSE4 loop.
  for (; length >= 4; length -= 4) {
    __m128 _t2 = _mm_loadu_ps(t);
    _t2 = _mm_max_ps(kZero, _mm_min_ps(_maxIdx, _t2));
    __m128 _w = _mm_round_ps(_t2, _MM_FROUND_FLOOR);

    idx1.v = _mm_cvtps_epi32(_w);
    idx2.v = _mm_cvtps_epi32(_mm_min_ps(_maxIdx, _mm_add_ps(_w, kOne)));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 4; ++k) {
      p1.f[k] = x[idx1.i[k]];
      p2.f[k] = x[idx2.i[k]];
    }

    _w = _mm_sub_ps(_t2, _w);
    _mm_storeu_ps(dst, _mm_add_ps(p1.v, _mm_mul_ps(_w, _mm_sub_ps(p2.v, p1.v))));

    dst += 4; t += 4;
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

void ArrayMathSSE4::sampleLinearRepeat_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (AM_UNLIKELY(xLength == 0)) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    ArrayMathSSE::fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32 xLengthF = static_cast<float32>(xLength);
  float32 xLengthFInv = 1.0f / xLengthF;
  __m128 _maxIdx = _mm_set1_ps(static_cast<float32>(maxIdx));
  __m128 _xLengthF = _mm_set1_ps(xLengthF);
  __m128 _xLengthFInv = _mm_set1_ps(xLengthFInv);
  const __m128 kOne = _mm_set1_ps(1.0f);

  union {
    __m128i v;
    int i[4];
  } idx1, idx2;

  union {
    __m128 v;
    float f[4];
  } p1, p2;

  // 1) Main SSE4 loop.
  for (; length >= 4; length -= 4) {
    __m128 _t2 = _mm_loadu_ps(t);
    _t2 = _mm_sub_ps(_t2, _mm_mul_ps(_mm_round_ps(_mm_mul_ps(_t2, _xLengthFInv), _MM_FROUND_FLOOR), _xLengthF));
    __m128 _w = _mm_round_ps(_t2, _MM_FROUND_FLOOR);

    idx1.v = _mm_cvtps_epi32(_w);
    idx2.v = _mm_cvtps_epi32(_mm_and_ps(_mm_cmpneq_ps(_w, _maxIdx), _mm_add_ps(_w, kOne)));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 4; ++k) {
      p1.f[k] = x[idx1.i[k]];
      p2.f[k] = x[idx2.i[k]];
    }

    _w = _mm_sub_ps(_t2, _w);
    _mm_storeu_ps(dst, _mm_add_ps(p1.v, _mm_mul_ps(_w, _mm_sub_ps(p2.v, p1.v))));

    dst += 4; t += 4;
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

#endif // AM_USE_X86 && AM_HAS_SSE4
