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

#include "x86/ArrayMathAVX.h"

#if defined(AM_USE_X86) && defined(AM_HAS_AVX)

#include <immintrin.h>

#include <cmath>

namespace arraymath {

//-----------------------------------------------------------------------------
// Template functions for common code patterns.
//-----------------------------------------------------------------------------

template <class OP>
void op_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  // 1) Align y to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(y) & 31) && length--) {
    *dst++ = OP::op(x, *y++);
  }
  const __m256 *_y = reinterpret_cast<const __m256*>(y);

  // Check alignment.
  bool dstAligned = (reinterpret_cast<size_t>(dst) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  __m256 _x = _mm256_set1_ps(x);
  if (dstAligned) {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_x, *_y++));
      dst += 8;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm256_storeu_ps(dst, OP::opAVX(_x, *_y++));
      dst += 8;
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
  // 1) Align x to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(x) & 31) && length--) {
    *dst++ = OP::op(*x++, *y++);
  }
  const __m256 *_x = reinterpret_cast<const __m256*>(x);

  // Check alignment.
  bool dstAligned = (reinterpret_cast<size_t>(dst) & 31) == 0;
  bool yAligned = (reinterpret_cast<size_t>(y) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  if (dstAligned && yAligned) {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(*_x++, _mm256_load_ps(y)));
      dst += 8;
      y += 8;
    }
  }
  else if (dstAligned) {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(*_x++, _mm256_loadu_ps(y)));
      dst += 8;
      y += 8;
    }
  }
  else if (yAligned) {
    for (; length >= 8; length -= 8) {
      _mm256_storeu_ps(dst, OP::opAVX(*_x++, _mm256_load_ps(y)));
      dst += 8;
      y += 8;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm256_storeu_ps(dst, OP::opAVX(*_x++, _mm256_loadu_ps(y)));
      dst += 8;
      y += 8;
    }
  }

  // 3) Tail loop.
  x = reinterpret_cast<const float32*>(_x);
  while (length--) {
    *dst++ = OP::op(*x++, *y++);
  }
}

template <class OP>
void op_f32_a(float32 *dst, const float32 *x, size_t length) {
  // 1) Align x to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(x) & 31) && length--) {
    *dst++ = OP::op(*x++);
  }

  // Check alignment.
  bool dstAligned = (reinterpret_cast<size_t>(dst) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  if (dstAligned) {
    for (; length >= 8; length -= 8) {
      _mm256_store_ps(dst, OP::opAVX(_mm256_load_ps(x)));
      dst += 8; x += 8;
    }
  }
  else {
    for (; length >= 8; length -= 8) {
      _mm256_storeu_ps(dst, OP::opAVX(_mm256_load_ps(x)));
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
    return _mm256_round_ps(a, _MM_FROUND_RINT);
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
  bool yAligned = (reinterpret_cast<size_t>(y) & 31) == 0;
  bool zAligned = (reinterpret_cast<size_t>(z) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  __m256 _x = _mm256_set1_ps(x);
  if (yAligned && zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_x, _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; y += 8; z += 8;
    }
  }
  else if (yAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_x, _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
      dst += 8; y += 8; z += 8;
    }
  }
  else if (zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_x, _mm256_loadu_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; y += 8; z += 8;
    }
  }
  else {
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
  bool xAligned = (reinterpret_cast<size_t>(x) & 31) == 0;
  bool yAligned = (reinterpret_cast<size_t>(y) & 31) == 0;
  bool zAligned = (reinterpret_cast<size_t>(z) & 31) == 0;

  // 2) Main AVX loop (handle different alignment cases).
  if (xAligned && yAligned && zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (xAligned && yAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (xAligned && zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_loadu_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (yAligned && zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (xAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_load_ps(x), _mm256_loadu_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (yAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_load_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_loadu_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else if (zAligned) {
    for (; length >= 8; length -= 8) {
      __m256 prod = _mm256_mul_ps(_mm256_loadu_ps(x), _mm256_loadu_ps(y));
      _mm256_store_ps(dst, _mm256_add_ps(prod, _mm256_load_ps(z)));
      dst += 8; x += 8; y += 8; z += 8;
    }
  }
  else {
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

void ArrayMathAVX::ceil_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CeilOP>(dst, x, length);
}

void ArrayMathAVX::floor_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FloorOP>(dst, x, length);
}

void ArrayMathAVX::round_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<RoundOP>(dst, x, length);
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_AVX
