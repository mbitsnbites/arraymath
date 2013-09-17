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

namespace {

#include "arm/neon_mathfun.h"

}

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

  // 2) Main NEON loop.
  float32x4_t _x = vdupq_n_f32(x);
  for (; length >= 4; length -= 4) {
    vst1q_f32(dst, OP::opNEON(_x, vld1q_f32(y)));
    dst += 4; y += 4;
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

  // 2) Main NEON loop.
  for (; length >= 4; length -= 4) {
    vst1q_f32(dst, OP::opNEON(vld1q_f32(x), vld1q_f32(y)));
    dst += 4; x += 4; y += 4;
  }

  // 3) Tail loop.
  while (length--) {
    *dst++ = OP::op(*x++, *y++);
  }
}

template <class OP>
void op_f32_a(float32 *dst, const float32 *x, size_t length) {
  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = OP::op(*x++);
  }

  // 2) Main NEON loop.
  for (; length >= 4; length -= 4) {
    vst1q_f32(dst, OP::opNEON(vld1q_f32(x)));
    dst += 4; x += 4;
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

void ArrayMathNEON::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  op_f32_sa<DivOP>(dst, x, y, length);
}

void ArrayMathNEON::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  op_f32_aa<DivOP>(dst, x, y, length);
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
  while ((reinterpret_cast<size_t>(x) & 15) && length--) {
    result = std::max(result, *x++);
  }

  // 2) Main NEON loop.
  float32x4_t _result = vdupq_n_f32(result);
  for (; length >= 4; length -= 4) {
    _result = vmaxq_f32(_result, vld1q_f32(x));
    x += 4;
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
  while ((reinterpret_cast<size_t>(x) & 15) && length--) {
    result = std::min(result, *x++);
  }

  // 2) Main NEON loop.
  float32x4_t _result = vdupq_n_f32(result);
  for (; length >= 4; length -= 4) {
    _result = vminq_f32(_result, vld1q_f32(x));
    x += 4;
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

void ArrayMathNEON::ramp_f32(float32 *dst, float32 first, float32 last, size_t length) {
  if (length == 0) {
    return;
  }
  if (length == 1) {
    *dst = first;
    return;
  }

  float32 step = (last - first) / static_cast<float32>(length - 1);
  float32 k = 0.0f;

  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst) & 15) && length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }

  // 2) Main NEON loop.
  static const float32x4_t kFour = vdupq_n_f32(4.0f);
  static const float32x4_t kIdxRamp = { 0.0f, 1.0f, 2.0f, 3.0f };
  float32x4_t _first = vdupq_n_f32(first);
  float32x4_t _step = vdupq_n_f32(step);
  float32x4_t _k = vaddq_f32(vdupq_n_f32(k), kIdxRamp);
  size_t mainLoopSize = (length / 4) * 4;
  for (; length >= 4; length -= 4) {
    vst1q_f32(dst, vaddq_f32(_first, vmulq_f32(_step, _k)));
    _k = vaddq_f32(_k, kFour);
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

#endif // AM_USE_ARM && AM_HAS_NEON
