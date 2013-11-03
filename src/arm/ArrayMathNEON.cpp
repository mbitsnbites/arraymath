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

#include "arm/ArrayMathNEON.h"

#if defined(AM_USE_ARM) && defined(AM_HAS_NEON)

#include <arm_neon.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace arraymath {

namespace {

#include "arm/neon_mathfun.h"


//-----------------------------------------------------------------------------
// Constants.
//-----------------------------------------------------------------------------

namespace constants {

const float32x4_t zero = { 0.0f, 0.0f, 0.0f, 0.0f };
const float32x4_t half = { 0.5f, 0.5f, 0.5f, 0.5f };
const float32x4_t one = { 1.0f, 1.0f, 1.0f, 1.0f };
const float32x4_t four = { 4.0f, 4.0f, 4.0f, 4.0f };
const float32x4_t round32 = { 12582912.0f, 12582912.0f, 12582912.0f, 12582912.0f };

} // namespace constants


//-----------------------------------------------------------------------------
// Aligned load/store operations.
//-----------------------------------------------------------------------------

AM_INLINE float32x4_t vld1qa_inc_f32(const float *&addr) {
  float32x4_t result;
#if defined(__GNUC__)
  __asm__ (
    "vld1.32 {%q[result]}, [%[addr]:128]!"
    : [result] "=w" (result), [addr] "+r" (addr)
    :
  );
#else
  result = vld1q_f32(addr);
  addr += 4;
#endif
  return result;
}

AM_INLINE void vst1qa_inc_f32(float *&addr, const float32x4_t x) {
#if defined(__GNUC__)
  __asm__ (
    "vst1.32 {%q[x]}, [%[addr]:128]!"
    : [addr] "+r" (addr)
    : [x] "w" (x)
  );
#else
  vst1q_f32(addr, x);
  addr += 4;
#endif
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

  // 2) Main NEON loop (handle different alignment cases).
  float32x4_t _x = vdupq_n_f32(x);
  if (aligned) {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, OP::opNEON(_x, vld1qa_inc_f32(y)));
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, OP::opNEON(_x, vld1q_f32(y)));
      y += 4;
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

  // 2) Main NEON loop.
  for (; length >= 4; length -= 4) {
    vst1qa_inc_f32(dst, OP::opNEON(vld1q_f32(x), vld1q_f32(y)));
    x += 4; y += 4;
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

  // 2) Main NEON loop (handle different alignment cases).
  if (aligned) {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, OP::opNEON(vld1qa_inc_f32(x)));
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, OP::opNEON(vld1q_f32(x)));
      x += 4;
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
  static float32x4_t opNEON(float32x4_t a, float32x4_t b) {
    return vaddq_f32(a, b);
  }
};

struct SubOP {
  static float32 op(float32 a, float32 b) {
    return a - b;
  }
  static float32x4_t opNEON(float32x4_t a, float32x4_t b) {
    return vsubq_f32(a, b);
  }
};

struct MulOP {
  static float32 op(float32 a, float32 b) {
    return a * b;
  }
  static float32x4_t opNEON(float32x4_t a, float32x4_t b) {
    return vmulq_f32(a, b);
  }
};

struct DivOP {
  static float32 op(float32 a, float32 b) {
    return a / b;
  }
  static float32x4_t opNEON(float32x4_t a, float32x4_t b) {
    // See http://rcl-rs-vvg.blogspot.se/2010_08_08_archive.html
    const float32x4_t q_inv0 = vrecpeq_f32(b);
    const float32x4_t q_step0 = vrecpsq_f32(q_inv0, b);
    const float32x4_t q_inv1 = vmulq_f32(q_step0, q_inv0);
    return vmulq_f32(a, q_inv1);
  }
};

struct SqrtOP {
  static float32 op(float32 a) {
    return std::sqrt(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    // See http://rcl-rs-vvg.blogspot.se/2010_08_08_archive.html
    const float32x4_t q_step_0 = vrsqrteq_f32(a);
    const float32x4_t q_step_parm0 = vmulq_f32(a, q_step_0);
    const float32x4_t q_step_result0 = vrsqrtsq_f32(q_step_parm0, q_step_0);
    const float32x4_t q_step_1 = vmulq_f32(q_step_0, q_step_result0);
    const float32x4_t q_step_parm1 = vmulq_f32(a, q_step_1);
    const float32x4_t q_step_result1 = vrsqrtsq_f32(q_step_parm1, q_step_1);
    const float32x4_t q_step_2 = vmulq_f32(q_step_1, q_step_result1);
    return vmulq_f32(a, q_step_2);
  }
};

struct AbsOP {
  static float32 op(float32 a) {
    return std::abs(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return vabsq_f32(a);
  }
};

struct CeilOP {
  static float32 op(float32 a) {
    return std::ceil(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    // This is essentially the same as RoundOP::opNEON(), but with two extras:
    // 1) We check if a is already an integer - in which case we keep it.
    // 2) Before rounding, we add 0.5 to a.
    float32x4_t rounded = vsubq_f32(vaddq_f32(a, constants::round32), constants::round32);
    uint32x4_t mask = vceqq_f32(a, rounded);
    float32x4_t ceiled = vsubq_f32(vaddq_f32(vaddq_f32(a, constants::half), constants::round32), constants::round32);
    return vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(a), mask), vbicq_u32(vreinterpretq_u32_f32(ceiled), mask)));
  }
};

struct FloorOP {
  static float32 op(float32 a) {
    return std::floor(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    // This is essentially the same as RoundOP::opNEON(), but with two extras:
    // 1) We check if a is already an integer - in which case we keep it.
    // 2) Before rounding, we subtract 0.5 from a.
    float32x4_t rounded = vsubq_f32(vaddq_f32(a, constants::round32), constants::round32);
    uint32x4_t mask = vceqq_f32(a, rounded);
    float32x4_t floored = vsubq_f32(vaddq_f32(vsubq_f32(a, constants::half), constants::round32), constants::round32);
    return vreinterpretq_f32_u32(vorrq_u32(vandq_u32(vreinterpretq_u32_f32(a), mask), vbicq_u32(vreinterpretq_u32_f32(floored), mask)));
  }
};

struct RoundOP {
  static float32 op(float32 a) {
    return std::floor(a + 0.5f);
  }
  static float32x4_t opNEON(float32x4_t a) {
    // This trick, (a + 1.5*2^23) - 1.5*2^23, flushes all the sub-integer bits,
    // using round to nearest mode (the only rounding mode supported by NEON).
    return vsubq_f32(vaddq_f32(a, constants::round32), constants::round32);
  }
};

struct SinOP {
  static float32 op(float32 a) {
    return std::sin(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return sin_ps(a);
  }
};

struct CosOP {
  static float32 op(float32 a) {
    return std::cos(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return cos_ps(a);
  }
};

struct ExpOP {
  static float32 op(float32 a) {
    return std::exp(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return exp_ps(a);
  }
};

struct LogOP {
  static float32 op(float32 a) {
    return std::log(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return log_ps(a);
  }
};

struct FractOP {
  static float32 op(float32 a) {
    return a - std::floor(a);
  }
  static float32x4_t opNEON(float32x4_t a) {
    return vsubq_f32(a, FloorOP::opNEON(a));
  }
};


} // anonymous namespace


//-----------------------------------------------------------------------------
// Class methods.
//-----------------------------------------------------------------------------

void ArrayMathNEON::add_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<AddOP>(dst, x, y, length);
}

void ArrayMathNEON::add_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<AddOP>(dst, x, y, length);
}

void ArrayMathNEON::sub_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<SubOP>(dst, x, y, length);
}

void ArrayMathNEON::sub_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<SubOP>(dst, x, y, length);
}

void ArrayMathNEON::mul_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<MulOP>(dst, x, y, length);
}

void ArrayMathNEON::mul_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<MulOP>(dst, x, y, length);
}

void ArrayMathNEON::mulCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main NEON loop.
  float32x4_t _xr = vdupq_n_f32(xReal);
  float32x4_t _xi = vdupq_n_f32(xImag);
  for (; length >= 12; length -= 12) {
    float32x4_t _yr1 = vld1q_f32(yReal);
    float32x4_t _yi1 = vld1q_f32(yImag);
    float32x4_t _yr2 = vld1q_f32(yReal + 4);
    float32x4_t _yi2 = vld1q_f32(yImag + 4);
    float32x4_t _yr3 = vld1q_f32(yReal + 8);
    float32x4_t _yi3 = vld1q_f32(yImag + 8);
    float32x4_t _xr_yr1 = vmulq_f32(_xr, _yr1);
    float32x4_t _xr_yi1 = vmulq_f32(_xr, _yi1);
    float32x4_t _xr_yr2 = vmulq_f32(_xr, _yr2);
    float32x4_t _xr_yi2 = vmulq_f32(_xr, _yi2);
    float32x4_t _xr_yr3 = vmulq_f32(_xr, _yr3);
    float32x4_t _xr_yi3 = vmulq_f32(_xr, _yi3);
    vst1q_f32(dstReal, vmlsq_f32(_xr_yr1, _xi, _yi1));
    vst1q_f32(dstImag, vmlaq_f32(_xr_yi1, _xi, _yr1));
    vst1q_f32(dstReal + 4, vmlsq_f32(_xr_yr2, _xi, _yi2));
    vst1q_f32(dstImag + 4, vmlaq_f32(_xr_yi2, _xi, _yr2));
    vst1q_f32(dstReal + 8, vmlsq_f32(_xr_yr3, _xi, _yi3));
    vst1q_f32(dstImag + 8, vmlaq_f32(_xr_yi3, _xi, _yr3));
    dstReal += 12; dstImag += 12; yReal += 12; yImag += 12;
  }
  for (; length >= 4; length -= 4) {
    float32x4_t _yr = vld1q_f32(yReal);
    float32x4_t _yi = vld1q_f32(yImag);
    float32x4_t _xr_yr = vmulq_f32(_xr, _yr);
    float32x4_t _xr_yi = vmulq_f32(_xr, _yi);
    vst1q_f32(dstReal, vmlsq_f32(_xr_yr, _xi, _yi));
    vst1q_f32(dstImag, vmlaq_f32(_xr_yi, _xi, _yr));
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

void ArrayMathNEON::mulCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  // 1) Main NEON loop.
  for (; length >= 12; length -= 12) {
      float32x4_t _xr1 = vld1q_f32(xReal);
      float32x4_t _xi1 = vld1q_f32(xImag);
      float32x4_t _yr1 = vld1q_f32(yReal);
      float32x4_t _yi1 = vld1q_f32(yImag);
      float32x4_t _xr2 = vld1q_f32(xReal + 4);
      float32x4_t _xi2 = vld1q_f32(xImag + 4);
      float32x4_t _yr2 = vld1q_f32(yReal + 4);
      float32x4_t _yi2 = vld1q_f32(yImag + 4);
      float32x4_t _xr3 = vld1q_f32(xReal + 8);
      float32x4_t _xi3 = vld1q_f32(xImag + 8);
      float32x4_t _yr3 = vld1q_f32(yReal + 8);
      float32x4_t _yi3 = vld1q_f32(yImag + 8);
      float32x4_t _xr_yr1 = vmulq_f32(_xr1, _yr1);
      float32x4_t _xr_yi1 = vmulq_f32(_xr1, _yi1);
      float32x4_t _xr_yr2 = vmulq_f32(_xr2, _yr2);
      float32x4_t _xr_yi2 = vmulq_f32(_xr2, _yi2);
      float32x4_t _xr_yr3 = vmulq_f32(_xr3, _yr3);
      float32x4_t _xr_yi3 = vmulq_f32(_xr3, _yi3);
      vst1q_f32(dstReal, vmlsq_f32(_xr_yr1, _xi1, _yi1));
      vst1q_f32(dstImag, vmlaq_f32(_xr_yi1, _xi1, _yr1));
      vst1q_f32(dstReal + 4, vmlsq_f32(_xr_yr2, _xi2, _yi2));
      vst1q_f32(dstImag + 4, vmlaq_f32(_xr_yi2, _xi2, _yr2));
      vst1q_f32(dstReal + 8, vmlsq_f32(_xr_yr3, _xi3, _yi3));
      vst1q_f32(dstImag + 8, vmlaq_f32(_xr_yi3, _xi3, _yr3));
      dstReal += 12; dstImag += 12; xReal += 12; xImag += 12; yReal += 12; yImag += 12;
  }
  for (; length >= 4; length -= 4) {
    float32x4_t _xr = vld1q_f32(xReal);
    float32x4_t _yr = vld1q_f32(yReal);
    float32x4_t _xi = vld1q_f32(xImag);
    float32x4_t _yi = vld1q_f32(yImag);
    float32x4_t _xr_yr = vmulq_f32(_xr, _yr);
    float32x4_t _xr_yi = vmulq_f32(_xr, _yi);
    vst1q_f32(dstReal, vmlsq_f32(_xr_yr, _xi, _yi));
    vst1q_f32(dstImag, vmlaq_f32(_xr_yi, _xi, _yr));
    dstReal += 4; dstImag += 4; xReal += 4; xImag += 4; yReal += 4; yImag += 4;
  }

  // 2) Tail loop.
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathNEON::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<DivOP>(dst, x, y, length);
}

void ArrayMathNEON::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<DivOP>(dst, x, y, length);
}

void ArrayMathNEON::madd_f32_saa(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = x * *y++ + *z++;
  }

  // 2) Main NEON loop.
  float32x4_t _x = vdupq_n_f32(x);
  for (; length >= 4; length -= 4) {
    float32x4_t _y = vld1q_f32(y);
    float32x4_t _z = vld1q_f32(z);
    vst1qa_inc_f32(dst, vmlaq_f32(_z, _x, _y));
    y += 4; z += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = x * *y++ + *z++;
  }
}

void ArrayMathNEON::madd_f32_aaa(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = *x++ * *y++ + *z++;
  }

  // 2) Main NEON loop.
  for (; length >= 4; length -= 4) {
    float32x4_t _x = vld1q_f32(x);
    float32x4_t _y = vld1q_f32(y);
    float32x4_t _z = vld1q_f32(z);
    vst1qa_inc_f32(dst, vmlaq_f32(_z, _x, _y));
    x += 4; y += 4; z += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = *x++ * *y++ + *z++;
  }
}

void ArrayMathNEON::abs_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<AbsOP>(dst, x, length);
}

void ArrayMathNEON::ceil_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CeilOP>(dst, x, length);
}

void ArrayMathNEON::floor_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FloorOP>(dst, x, length);
}

void ArrayMathNEON::round_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<RoundOP>(dst, x, length);
}

void ArrayMathNEON::sin_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<SinOP>(dst, x, length);
}

void ArrayMathNEON::cos_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<CosOP>(dst, x, length);
}

void ArrayMathNEON::exp_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<ExpOP>(dst, x, length);
}

void ArrayMathNEON::log_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<LogOP>(dst, x, length);
}

void ArrayMathNEON::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<SqrtOP>(dst, x, length);
}

float32 ArrayMathNEON::max_f32(const float32 *x, size_t length) {
  float32 result = -std::numeric_limits<float>::infinity();

  // 1) Align x to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(x) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    result = std::max(result, *x++);
  }

  // 2) Main NEON loop.
  float32x4_t _result = vdupq_n_f32(result);
  for (; length >= 4; length -= 4) {
    _result = vmaxq_f32(_result, vld1qa_inc_f32(x));
  }

  // 3) Horizontal max of SIMD register (there must be a better way?).
  union {
    float32x4_t v;
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

float32 ArrayMathNEON::min_f32(const float32 *x, size_t length) {
  float32 result = std::numeric_limits<float>::infinity();

  // 1) Align x to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(x) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    result = std::min(result, *x++);
  }

  // 2) Main NEON loop.
  float32x4_t _result = vdupq_n_f32(result);
  for (; length >= 4; length -= 4) {
    _result = vminq_f32(_result, vld1qa_inc_f32(x));
  }

  // 3) Horizontal min of SIMD register (there must be a better way?).
  union {
    float32x4_t v;
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

void ArrayMathNEON::clamp_f32(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
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

  // 2) Main NEON loop (handle different alignment cases).
  float32x4_t _xMin = vdupq_n_f32(xMin);
  float32x4_t _xMax = vdupq_n_f32(xMax);
  if (aligned) {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, vmaxq_f32(_xMin, vminq_f32(_xMax, vld1qa_inc_f32(x))));
    }
  }
  else {
    for (; length >= 4; length -= 4) {
      vst1qa_inc_f32(dst, vmaxq_f32(_xMin, vminq_f32(_xMax, vld1q_f32(x))));
      x += 4;
    }
  }

  // 3) Tail loop.
  while (length--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }
}

void ArrayMathNEON::fract_f32(float32 *dst, const float32 *x, size_t length) {
  op_f32_a<FractOP>(dst, x, length);
}

void ArrayMathNEON::fill_f32(float32 *dst, float32 value, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  size_t numUnaligned = (reinterpret_cast<size_t>(dst) & 15) >> 2;
  numUnaligned = std::min(numUnaligned ? 4 - numUnaligned : 0, length);
  length -= numUnaligned;
  while (numUnaligned--) {
    *dst++ = value;
  }

  // 2) Main NEON loop.
  float32x4_t _value = vdupq_n_f32(value);
  for (; length >= 4; length -= 4) {
    vst1qa_inc_f32(dst, _value);
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = value;
  }
}

void ArrayMathNEON::ramp_f32(float32 *dst, float32 first, float32 last, size_t length) {
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

  // 2) Main NEON loop.
  static const float32x4_t kIdxRamp = { 0.0f, 1.0f, 2.0f, 3.0f };
  float32x4_t _first = vdupq_n_f32(first);
  float32x4_t _step = vdupq_n_f32(step);
  float32x4_t _k = vaddq_f32(vdupq_n_f32(k), kIdxRamp);
  size_t mainLoopSize = (length / 4) * 4;
  for (; length >= 4; length -= 4) {
    vst1qa_inc_f32(dst, vaddq_f32(_first, vmulq_f32(_step, _k)));
    _k = vaddq_f32(_k, constants::four);
  }
  k += mainLoopSize;

  // 3) Tail loop.
  while (length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }
}

void ArrayMathNEON::sampleLinear_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (AM_UNLIKELY(xLength == 0)) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32x4_t _maxIdx = vdupq_n_f32(static_cast<float32>(maxIdx));

  union {
    uint32x4_t v;
    unsigned int u[4];
  } idx1, idx2;

  union {
    float32x4_t v;
    int i[4];
  } p1, p2;

  // 1) Main NEON loop.
  for (; length >= 4; length -= 4) {
    float32x4_t _t2 = vld1q_f32(t);
    _t2 = vmaxq_f32(constants::zero, vminq_f32(_maxIdx, _t2));
    float32x4_t _w = FloorOP::opNEON(_t2);

    idx1.v = vcvtq_u32_f32(_w);
    idx2.v = vcvtq_u32_f32(vminq_f32(_maxIdx, vaddq_f32(_w, constants::one)));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 4; ++k) {
      p1.i[k] = reinterpret_cast<const int*>(x)[idx1.u[k]];
      p2.i[k] = reinterpret_cast<const int*>(x)[idx2.u[k]];
    }

    _w = vsubq_f32(_t2, _w);
    vst1q_f32(dst, vmlaq_f32(p1.v, _w, vsubq_f32(p2.v, p1.v)));

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

void ArrayMathNEON::sampleLinearRepeat_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (AM_UNLIKELY(xLength == 0)) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32 xLengthF = static_cast<float32>(xLength);
  float32 xLengthFInv = 1.0f / xLengthF;
  float32x4_t _maxIdx = vdupq_n_f32(static_cast<float32>(maxIdx));
  float32x4_t _xLengthF = vdupq_n_f32(xLengthF);
  float32x4_t _xLengthFInv = vdupq_n_f32(xLengthFInv);

  union {
    uint32x4_t v;
    unsigned int u[4];
  } idx1, idx2;

  union {
    float32x4_t v;
    int i[4];
  } p1, p2;

  // 1) Main NEON loop.
  for (; length >= 4; length -= 4) {
    float32x4_t _t2 = vld1q_f32(t);
    _t2 = vsubq_f32(_t2, vmulq_f32(FloorOP::opNEON(vmulq_f32(_t2, _xLengthFInv)), _xLengthF));
    float32x4_t _w = FloorOP::opNEON(_t2);

    idx1.v = vcvtq_u32_f32(_w);
    idx2.v = vbicq_u32(vcvtq_u32_f32(vaddq_f32(_w, constants::one)), vceqq_f32(_w, _maxIdx));

    // TODO(m): Can we do this in a better way?
    for (int k = 0; k < 4; ++k) {
      p1.i[k] = reinterpret_cast<const int*>(x)[idx1.u[k]];
      p2.i[k] = reinterpret_cast<const int*>(x)[idx2.u[k]];
    }

    _w = vsubq_f32(_t2, _w);
    vst1q_f32(dst, vmlaq_f32(p1.v, _w, vsubq_f32(p2.v, p1.v)));

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

#endif // AM_USE_ARM && AM_HAS_NEON
