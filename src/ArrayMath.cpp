// -*- Mode: c++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//------------------------------------------------------------------------------
// ArrayMath - an array math library
//------------------------------------------------------------------------------
// Copyright(c) 2013 Marcus Geelnard
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

#include "Architecture.h"
#include "ArrayMath.h"
#include "ArrayMathGeneric.h"
#include "ArrayMathSSE.h"
#include "CPUFeatureDetector.h"
#include "RandomGeneric.h"

namespace arraymath {

ArrayMath::ArrayMath() {
  // Set up pointers to generic routines.
  p_add_f32_sa = ArrayMathGeneric::add_f32_sa;
  p_add_f32_aa = ArrayMathGeneric::add_f32_aa;
  p_sub_f32_sa = ArrayMathGeneric::sub_f32_sa;
  p_sub_f32_aa = ArrayMathGeneric::sub_f32_aa;
  p_mul_f32_sa = ArrayMathGeneric::mul_f32_sa;
  p_mul_f32_aa = ArrayMathGeneric::mul_f32_aa;
  p_mulCplx_f32_sa = ArrayMathGeneric::mulCplx_f32_sa;
  p_mulCplx_f32_aa = ArrayMathGeneric::mulCplx_f32_aa;
  p_div_f32_sa = ArrayMathGeneric::div_f32_sa;
  p_div_f32_aa = ArrayMathGeneric::div_f32_aa;
  p_divCplx_f32_sa = ArrayMathGeneric::divCplx_f32_sa;
  p_divCplx_f32_aa = ArrayMathGeneric::divCplx_f32_aa;
  p_madd_f32_saa = ArrayMathGeneric::madd_f32_saa;
  p_madd_f32_aaa = ArrayMathGeneric::madd_f32_aaa;
  p_abs_f32 = ArrayMathGeneric::abs_f32;
  p_absCplx_f32 = ArrayMathGeneric::absCplx_f32;
  p_acos_f32 = ArrayMathGeneric::acos_f32;
  p_asin_f32 = ArrayMathGeneric::asin_f32;
  p_atan_f32 = ArrayMathGeneric::atan_f32;
  p_atan2_f32 = ArrayMathGeneric::atan2_f32;
  p_ceil_f32 = ArrayMathGeneric::ceil_f32;
  p_cos_f32 = ArrayMathGeneric::cos_f32;
  p_exp_f32 = ArrayMathGeneric::exp_f32;
  p_floor_f32 = ArrayMathGeneric::floor_f32;
  p_log_f32 = ArrayMathGeneric::log_f32;
  p_max_f32 = ArrayMathGeneric::max_f32;
  p_min_f32 = ArrayMathGeneric::min_f32;
  p_pow_f32_as = ArrayMathGeneric::pow_f32_as;
  p_pow_f32_aa = ArrayMathGeneric::pow_f32_aa;
  p_round_f32 = ArrayMathGeneric::round_f32;
  p_sin_f32 = ArrayMathGeneric::sin_f32;
  p_sqrt_f32 = ArrayMathGeneric::sqrt_f32;
  p_tan_f32 = ArrayMathGeneric::tan_f32;
  p_clamp_f32 = ArrayMathGeneric::clamp_f32;
  p_fract_f32 = ArrayMathGeneric::fract_f32;
  p_ramp_f32 = ArrayMathGeneric::ramp_f32;
  p_sign_f32 = ArrayMathGeneric::sign_f32;
  p_sum_f32 = ArrayMathGeneric::sum_f32;
  p_sampleLinear_f32 = ArrayMathGeneric::sampleLinear_f32;
  p_sampleLinearRepeat_f32 = ArrayMathGeneric::sampleLinearRepeat_f32;
  p_sampleCubic_f32 = ArrayMathGeneric::sampleCubic_f32;
  p_sampleCubicRepeat_f32 = ArrayMathGeneric::sampleCubicRepeat_f32;

  // Set up default random number generator.
  m_random = new RandomGeneric();

#ifdef AM_USE_X86
  CPUFeatureDetector cpu;

#ifdef AM_HAS_SSE
  if (cpu.hasSSE()) {
    // Override generic routines with x86 SSE optimized versions.
    p_add_f32_sa = ArrayMathSSE::add_f32_sa;
    p_add_f32_aa = ArrayMathSSE::add_f32_aa;
    p_sub_f32_sa = ArrayMathSSE::sub_f32_sa;
    p_sub_f32_aa = ArrayMathSSE::sub_f32_aa;
    p_mul_f32_sa = ArrayMathSSE::mul_f32_sa;
    p_mul_f32_aa = ArrayMathSSE::mul_f32_aa;
  }
#endif // AM_HAS_SSE

#endif // AM_USE_X86

#ifdef AM_USE_ARM
  CPUFeatureDetector cpu;

#ifdef AM_HAS_NEON
  if (cpu.hasNEON()) {
    // TODO(m): Override generic routines with ARM NEON optimized versions.
    if (cpu.hasNEON_FMA()) {
      // TODO(m): madd() should probably be implemented using vmla.
    }
  }
#endif // AM_HAS_NEON

#endif // AM_USE_ARM
}

ArrayMath::~ArrayMath() {
  delete m_random;
}

}  // namespace arraymath
