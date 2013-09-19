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

#include "ArrayMath.h"

#include "arm/ArrayMathNEON.h"
#include "common/Architecture.h"
#include "common/CPUFeatureDetector.h"
#include "generic/ArrayMathGeneric.h"
#include "generic/RandomGeneric.h"
#include "x86/ArrayMathAVX.h"
#include "x86/ArrayMathSSE.h"
#include "x86/ArrayMathSSE2.h"
#include "x86/RandomSSE2.h"

namespace arraymath {

ArrayMath::ArrayMath() {
  CPUFeatureDetector cpu;

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

#if defined(AM_USE_X86) && defined(AM_HAS_SSE)
  if (cpu.hasSSE()) {
    // Override generic routines with x86 SSE optimized versions.
    p_add_f32_sa = ArrayMathSSE::add_f32_sa;
    p_add_f32_aa = ArrayMathSSE::add_f32_aa;
    p_sub_f32_sa = ArrayMathSSE::sub_f32_sa;
    p_sub_f32_aa = ArrayMathSSE::sub_f32_aa;
    p_mul_f32_sa = ArrayMathSSE::mul_f32_sa;
    p_mul_f32_aa = ArrayMathSSE::mul_f32_aa;
    p_mulCplx_f32_sa = ArrayMathSSE::mulCplx_f32_sa;
    p_mulCplx_f32_aa = ArrayMathSSE::mulCplx_f32_aa;
    p_max_f32 = ArrayMathSSE::max_f32;
    p_min_f32 = ArrayMathSSE::min_f32;
    p_sqrt_f32 = ArrayMathSSE::sqrt_f32;
    p_sin_f32 = ArrayMathSSE::sin_f32;
    p_cos_f32 = ArrayMathSSE::cos_f32;
    p_exp_f32 = ArrayMathSSE::exp_f32;
    p_log_f32 = ArrayMathSSE::log_f32;
    p_ramp_f32 = ArrayMathSSE::ramp_f32;
  }
#endif // AM_USE_X86 && AM_HAS_SSE

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)
  if (cpu.hasSSE2()) {
    // Override generic routines with x86 SSE2 optimized versions.
    p_div_f32_sa = ArrayMathSSE2::div_f32_sa;
    p_div_f32_aa = ArrayMathSSE2::div_f32_aa;
  }
#endif // AM_USE_X86 && AM_HAS_SSE2

#if defined(AM_USE_X86) && defined(AM_HAS_AVX)
  if (cpu.hasAVX()) {
    // Override generic routines with x86 SSE optimized versions.
    p_add_f32_sa = ArrayMathAVX::add_f32_sa;
    p_add_f32_aa = ArrayMathAVX::add_f32_aa;
    p_sub_f32_sa = ArrayMathAVX::sub_f32_sa;
    p_sub_f32_aa = ArrayMathAVX::sub_f32_aa;
    p_mul_f32_sa = ArrayMathAVX::mul_f32_sa;
    p_mul_f32_aa = ArrayMathAVX::mul_f32_aa;
    p_mulCplx_f32_sa = ArrayMathAVX::mulCplx_f32_sa;
    p_mulCplx_f32_aa = ArrayMathAVX::mulCplx_f32_aa;
    p_div_f32_sa = ArrayMathAVX::div_f32_sa;
    p_div_f32_aa = ArrayMathAVX::div_f32_aa;
    p_madd_f32_saa = ArrayMathAVX::madd_f32_saa;
    p_madd_f32_aaa = ArrayMathAVX::madd_f32_aaa;
    p_sqrt_f32 = ArrayMathAVX::sqrt_f32;
    p_ceil_f32 = ArrayMathAVX::ceil_f32;
    p_floor_f32 = ArrayMathAVX::floor_f32;
    p_max_f32 = ArrayMathAVX::max_f32;
    p_min_f32 = ArrayMathAVX::min_f32;
    p_round_f32 = ArrayMathAVX::round_f32;
    p_clamp_f32 = ArrayMathAVX::clamp_f32;
    p_fract_f32 = ArrayMathAVX::fract_f32;
    p_ramp_f32 = ArrayMathAVX::ramp_f32;
    p_sign_f32 = ArrayMathAVX::sign_f32;
    p_sampleLinear_f32 = ArrayMathAVX::sampleLinear_f32;
  }
#endif // AM_USE_X86 && AM_HAS_AVX

#if defined(AM_USE_ARM) && defined(AM_HAS_NEON)
  if (cpu.hasNEON()) {
    p_add_f32_sa = ArrayMathNEON::add_f32_sa;
    p_add_f32_aa = ArrayMathNEON::add_f32_aa;
    p_sub_f32_sa = ArrayMathNEON::sub_f32_sa;
    p_sub_f32_aa = ArrayMathNEON::sub_f32_aa;
    p_mul_f32_sa = ArrayMathNEON::mul_f32_sa;
    p_mul_f32_aa = ArrayMathNEON::mul_f32_aa;
    p_mulCplx_f32_sa = ArrayMathNEON::mulCplx_f32_sa;
    p_mulCplx_f32_aa = ArrayMathNEON::mulCplx_f32_aa;
    p_div_f32_sa = ArrayMathNEON::div_f32_sa;
    p_div_f32_aa = ArrayMathNEON::div_f32_aa;
    p_sin_f32 = ArrayMathNEON::sin_f32;
    p_cos_f32 = ArrayMathNEON::cos_f32;
    p_exp_f32 = ArrayMathNEON::exp_f32;
    p_log_f32 = ArrayMathNEON::log_f32;
    p_sqrt_f32 = ArrayMathNEON::sqrt_f32;
    p_max_f32 = ArrayMathNEON::max_f32;
    p_min_f32 = ArrayMathNEON::min_f32;
    p_ramp_f32 = ArrayMathNEON::ramp_f32;
    if (cpu.hasNEON_FMA()) {
      // TODO(m): madd() should probably be implemented using vmla.
    }
  }
#endif // AM_USE_ARM && AM_HAS_NEON

  // Set up random number generator.
  m_random_f32 = NULL;

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)
  if (cpu.hasSSE2()) {
    // Use an SSE2 optimized random number generator.
    m_random_f32 = RandomSSE2Factory::create();
  }
#endif // AM_USE_X86 && AM_HAS_SSE2

  // Fall back to generic random number generator if necessary.
  if (!m_random_f32) {
    m_random_f32 = new RandomGeneric();
  }
}

ArrayMath::~ArrayMath() {
  delete m_random_f32;
}

}  // namespace arraymath
