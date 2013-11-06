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

#include "arm/ArrayMathARM.h"
#include "arm/ArrayMathNEON.h"
#include "arm/RandomNEON.h"
#include "common/Architecture.h"
#include "common/CPUFeatureDetector.h"
#include "generic/ArrayMathGeneric.h"
#include "generic/RandomGeneric.h"
#include "x86/ArrayMathAVX.h"
#include "x86/ArrayMathSSE.h"
#include "x86/ArrayMathSSE4.h"
#include "x86/RandomSSE2.h"

namespace arraymath {

ArrayMath::ArrayMath() {
  CPUFeatureDetector cpu;

  // Set up pointers to generic routines.
  m_add_f32_sa = ArrayMathGeneric::add_f32_sa;
  m_add_f32_aa = ArrayMathGeneric::add_f32_aa;
  m_sub_f32_sa = ArrayMathGeneric::sub_f32_sa;
  m_sub_f32_aa = ArrayMathGeneric::sub_f32_aa;
  m_mul_f32_sa = ArrayMathGeneric::mul_f32_sa;
  m_mul_f32_aa = ArrayMathGeneric::mul_f32_aa;
  m_mulCplx_f32_sa = ArrayMathGeneric::mulCplx_f32_sa;
  m_mulCplx_f32_aa = ArrayMathGeneric::mulCplx_f32_aa;
  m_div_f32_sa = ArrayMathGeneric::div_f32_sa;
  m_div_f32_aa = ArrayMathGeneric::div_f32_aa;
  m_divCplx_f32_sa = ArrayMathGeneric::divCplx_f32_sa;
  m_divCplx_f32_aa = ArrayMathGeneric::divCplx_f32_aa;
  m_madd_f32_saa = ArrayMathGeneric::madd_f32_saa;
  m_madd_f32_aaa = ArrayMathGeneric::madd_f32_aaa;
  m_abs_f32 = ArrayMathGeneric::abs_f32;
  m_absCplx_f32 = ArrayMathGeneric::absCplx_f32;
  m_acos_f32 = ArrayMathGeneric::acos_f32;
  m_asin_f32 = ArrayMathGeneric::asin_f32;
  m_atan_f32 = ArrayMathGeneric::atan_f32;
  m_atan2_f32 = ArrayMathGeneric::atan2_f32;
  m_ceil_f32 = ArrayMathGeneric::ceil_f32;
  m_cos_f32 = ArrayMathGeneric::cos_f32;
  m_exp_f32 = ArrayMathGeneric::exp_f32;
  m_floor_f32 = ArrayMathGeneric::floor_f32;
  m_log_f32 = ArrayMathGeneric::log_f32;
  m_max_f32 = ArrayMathGeneric::max_f32;
  m_min_f32 = ArrayMathGeneric::min_f32;
  m_pow_f32_as = ArrayMathGeneric::pow_f32_as;
  m_pow_f32_aa = ArrayMathGeneric::pow_f32_aa;
  m_round_f32 = ArrayMathGeneric::round_f32;
  m_sin_f32 = ArrayMathGeneric::sin_f32;
  m_sqrt_f32 = ArrayMathGeneric::sqrt_f32;
  m_tan_f32 = ArrayMathGeneric::tan_f32;
  m_clamp_f32 = ArrayMathGeneric::clamp_f32;
  m_fract_f32 = ArrayMathGeneric::fract_f32;
  m_fill_f32 = ArrayMathGeneric::fill_f32;
  m_ramp_f32 = ArrayMathGeneric::ramp_f32;
  m_sign_f32 = ArrayMathGeneric::sign_f32;
  m_sum_f32 = ArrayMathGeneric::sum_f32;
  m_sampleLinear_f32 = ArrayMathGeneric::sampleLinear_f32;
  m_sampleLinearRepeat_f32 = ArrayMathGeneric::sampleLinearRepeat_f32;
  m_sampleCubic_f32 = ArrayMathGeneric::sampleCubic_f32;
  m_sampleCubicRepeat_f32 = ArrayMathGeneric::sampleCubicRepeat_f32;

#if defined(AM_USE_X86) && defined(AM_HAS_SSE)
  if (cpu.hasSSE()) {
    // Override generic routines with x86 SSE optimized versions.
    m_add_f32_sa = ArrayMathSSE::add_f32_sa;
    m_add_f32_aa = ArrayMathSSE::add_f32_aa;
    m_sub_f32_sa = ArrayMathSSE::sub_f32_sa;
    m_sub_f32_aa = ArrayMathSSE::sub_f32_aa;
    m_mul_f32_sa = ArrayMathSSE::mul_f32_sa;
    m_mul_f32_aa = ArrayMathSSE::mul_f32_aa;
    m_mulCplx_f32_sa = ArrayMathSSE::mulCplx_f32_sa;
    m_mulCplx_f32_aa = ArrayMathSSE::mulCplx_f32_aa;
    m_div_f32_sa = ArrayMathSSE::div_f32_sa;
    m_div_f32_aa = ArrayMathSSE::div_f32_aa;
    m_divCplx_f32_sa = ArrayMathSSE::divCplx_f32_sa;
    m_divCplx_f32_aa = ArrayMathSSE::divCplx_f32_aa;
    m_max_f32 = ArrayMathSSE::max_f32;
    m_min_f32 = ArrayMathSSE::min_f32;
    m_sqrt_f32 = ArrayMathSSE::sqrt_f32;
    m_abs_f32 = ArrayMathSSE::abs_f32;
    m_absCplx_f32 = ArrayMathSSE::absCplx_f32;
    m_sin_f32 = ArrayMathSSE::sin_f32;
    m_cos_f32 = ArrayMathSSE::cos_f32;
    m_tan_f32 = ArrayMathSSE::tan_f32;
    m_exp_f32 = ArrayMathSSE::exp_f32;
    m_log_f32 = ArrayMathSSE::log_f32;
    m_pow_f32_as = ArrayMathSSE::pow_f32_as;
    m_pow_f32_aa = ArrayMathSSE::pow_f32_aa;
    m_clamp_f32 = ArrayMathSSE::clamp_f32;
    m_fill_f32 = ArrayMathSSE::fill_f32;
    m_ramp_f32 = ArrayMathSSE::ramp_f32;
  }
#endif // AM_USE_X86 && AM_HAS_SSE

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)
  if (cpu.hasSSE2()) {
    // TODO(m): Override generic routines with x86 SSE2 optimized versions.
  }
#endif // AM_USE_X86 && AM_HAS_SSE2

#if defined(AM_USE_X86) && defined(AM_HAS_SSE4)
  if (cpu.hasSSE4()) {
    m_ceil_f32 = ArrayMathSSE4::ceil_f32;
    m_floor_f32 = ArrayMathSSE4::floor_f32;
    m_round_f32 = ArrayMathSSE4::round_f32;
    m_fract_f32 = ArrayMathSSE4::fract_f32;
    m_sampleLinear_f32 = ArrayMathSSE4::sampleLinear_f32;
    m_sampleLinearRepeat_f32 = ArrayMathSSE4::sampleLinearRepeat_f32;
  }
#endif // AM_USE_X86 && AM_HAS_SSE4

#if defined(AM_USE_X86) && defined(AM_HAS_AVX)
  if (cpu.hasAVX()) {
    // Override generic routines with x86 SSE optimized versions.
    m_add_f32_sa = ArrayMathAVX::add_f32_sa;
    m_add_f32_aa = ArrayMathAVX::add_f32_aa;
    m_sub_f32_sa = ArrayMathAVX::sub_f32_sa;
    m_sub_f32_aa = ArrayMathAVX::sub_f32_aa;
    m_mul_f32_sa = ArrayMathAVX::mul_f32_sa;
    m_mul_f32_aa = ArrayMathAVX::mul_f32_aa;
    m_mulCplx_f32_sa = ArrayMathAVX::mulCplx_f32_sa;
    m_mulCplx_f32_aa = ArrayMathAVX::mulCplx_f32_aa;
    m_div_f32_sa = ArrayMathAVX::div_f32_sa;
    m_div_f32_aa = ArrayMathAVX::div_f32_aa;
    m_divCplx_f32_sa = ArrayMathAVX::divCplx_f32_sa;
    m_divCplx_f32_aa = ArrayMathAVX::divCplx_f32_aa;
    m_madd_f32_saa = ArrayMathAVX::madd_f32_saa;
    m_madd_f32_aaa = ArrayMathAVX::madd_f32_aaa;
    m_abs_f32 = ArrayMathAVX::abs_f32;
    m_absCplx_f32 = ArrayMathAVX::absCplx_f32;
    m_sqrt_f32 = ArrayMathAVX::sqrt_f32;
    m_ceil_f32 = ArrayMathAVX::ceil_f32;
    m_floor_f32 = ArrayMathAVX::floor_f32;
    m_max_f32 = ArrayMathAVX::max_f32;
    m_min_f32 = ArrayMathAVX::min_f32;
    m_round_f32 = ArrayMathAVX::round_f32;
    m_clamp_f32 = ArrayMathAVX::clamp_f32;
    m_fract_f32 = ArrayMathAVX::fract_f32;
    m_fill_f32 = ArrayMathAVX::fill_f32;
    m_ramp_f32 = ArrayMathAVX::ramp_f32;
    m_sign_f32 = ArrayMathAVX::sign_f32;
    m_sampleLinear_f32 = ArrayMathAVX::sampleLinear_f32;
    m_sampleLinearRepeat_f32 = ArrayMathAVX::sampleLinearRepeat_f32;
  }
#endif // AM_USE_X86 && AM_HAS_AVX

#if defined(AM_USE_ARM)
  m_abs_f32 = ArrayMathARM::abs_f32;
  m_fill_f32 = ArrayMathARM::fill_f32;
  m_sign_f32 = ArrayMathARM::sign_f32;
#endif // AM_USE_ARM

#if defined(AM_USE_ARM) && defined(AM_HAS_NEON)
  if (cpu.hasNEON()) {
    m_add_f32_sa = ArrayMathNEON::add_f32_sa;
    m_add_f32_aa = ArrayMathNEON::add_f32_aa;
    m_sub_f32_sa = ArrayMathNEON::sub_f32_sa;
    m_sub_f32_aa = ArrayMathNEON::sub_f32_aa;
    m_mul_f32_sa = ArrayMathNEON::mul_f32_sa;
    m_mul_f32_aa = ArrayMathNEON::mul_f32_aa;
    m_mulCplx_f32_sa = ArrayMathNEON::mulCplx_f32_sa;
    m_mulCplx_f32_aa = ArrayMathNEON::mulCplx_f32_aa;
    m_div_f32_sa = ArrayMathNEON::div_f32_sa;
    m_div_f32_aa = ArrayMathNEON::div_f32_aa;
    m_madd_f32_saa = ArrayMathNEON::madd_f32_saa;
    m_madd_f32_aaa = ArrayMathNEON::madd_f32_aaa;
    m_abs_f32 = ArrayMathNEON::abs_f32;
    m_sin_f32 = ArrayMathNEON::sin_f32;
    m_cos_f32 = ArrayMathNEON::cos_f32;
    m_exp_f32 = ArrayMathNEON::exp_f32;
    m_log_f32 = ArrayMathNEON::log_f32;
    m_sqrt_f32 = ArrayMathNEON::sqrt_f32;
    m_ceil_f32 = ArrayMathNEON::ceil_f32;
    m_floor_f32 = ArrayMathNEON::floor_f32;
    m_max_f32 = ArrayMathNEON::max_f32;
    m_min_f32 = ArrayMathNEON::min_f32;
    m_round_f32 = ArrayMathNEON::round_f32;
    m_clamp_f32 = ArrayMathNEON::clamp_f32;
    m_fract_f32 = ArrayMathNEON::fract_f32;
    m_fill_f32 = ArrayMathNEON::fill_f32;
    m_ramp_f32 = ArrayMathNEON::ramp_f32;
    m_sampleLinear_f32 = ArrayMathNEON::sampleLinear_f32;
    m_sampleLinearRepeat_f32 = ArrayMathNEON::sampleLinearRepeat_f32;
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

#if defined(AM_USE_ARM) && defined(AM_HAS_NEON)
  if (cpu.hasNEON()) {
    // Use a NEON optimized random number generator.
    m_random_f32 = RandomNEONFactory::create();
  }
#endif // AM_USE_ARM && AM_HAS_NEON

  // Fall back to generic random number generator if necessary.
  if (!m_random_f32) {
    m_random_f32 = new RandomGeneric();
  }
}

ArrayMath::~ArrayMath() {
  delete m_random_f32;
}

}  // namespace arraymath
