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
// TODO: We get a significant performance hit for the run in/out loops if this
// unit is compiled for x87 or SSE/SSE2 standard math (which is usually the
// default setting), which is painful for short arrays. This is due to the
// expensive switch between SSE & AVX. We should implement pure AVX solutions
// wherever possible!
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

void ArrayMathAVX::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 32-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 31) && length--) {
    *dst++ = std::sqrt(*x++);
  }

  // Check alignment.
  bool aligned = (reinterpret_cast<size_t>(x) & 31) == 0;

  // 2) Main SSE loop (handle different alignment cases).
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

void ArrayMathAVX::round_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<RoundOP>(dst, x, length);
}

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_AVX
