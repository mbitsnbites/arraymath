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

#ifndef _ARRAYMATH_ARRAYMATH_H
#define _ARRAYMATH_ARRAYMATH_H

#include "Random.h"
#include "Types.h"

namespace arraymath {

/// \brief Array math class.
///
/// ### Using the ArrayMath interface
///
/// \code
///    #include <ArrayMath.h>
///
///    void myFunction {
///      // Initialize the ArrayMath object. Note: This object may contain state (for
///      // instance for the random number generator), so it is *not* thread safe.
///      arraymath::ArrayMath math;
///
///      // Call array methods.
///      const unsigned len = 128;
///      float a[len], b[len];
///      math.ramp(a, 0.0f, 100.0f, len);    // a = [0.0 .. 100.0]
///      math.random(b, -1.0f, 1.0f, len);   // b = random, [-1.0, 1.0)
///      math.madd(a, 0.5f, a, b, len);      // a = 0.5 * a + b
///    }
/// \endcode
class ArrayMath {
 public:
  ArrayMath();
  ~ArrayMath();

  /// Compute x + y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source scalar.
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void add(float32 *dst, float32 x, const float32 *y, size_t length) {
    p_add_f32_sa(dst, x, y, length);
  }

  /// Compute x + y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void add(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_add_f32_aa(dst, x, y, length);
  }

  /// Compute x - y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source scalar.
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void sub(float32 *dst, float32 x, const float32 *y, size_t length) {
    p_sub_f32_sa(dst, x, y, length);
  }

  /// Compute x - y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void sub(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_sub_f32_aa(dst, x, y, length);
  }

  /// Compute x * y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source scalar.
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void mul(float32 *dst, float32 x, const float32 *y, size_t length) {
    p_mul_f32_sa(dst, x, y, length);
  }

  /// Compute x * y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void mul(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_mul_f32_aa(dst, x, y, length);
  }

  AM_INLINE void mulCplx(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_mulCplx_f32_sa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  AM_INLINE void mulCplx(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_mulCplx_f32_aa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  /// Compute x / y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source scalar.
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void div(float32 *dst, float32 x, const float32 *y, size_t length) {
    p_div_f32_sa(dst, x, y, length);
  }

  /// Compute x / y and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void div(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_div_f32_aa(dst, x, y, length);
  }

  AM_INLINE void divCplx(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_divCplx_f32_sa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  AM_INLINE void divCplx(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_divCplx_f32_aa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  /// Compute x * y + z and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source scalar.
  /// \param[in] y Source array (\c length elements).
  /// \param[in] z Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void madd(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length) {
    p_madd_f32_saa(dst, x, y, z, length);
  }

  /// Compute x * y + z and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] z Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void madd(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length) {
    p_madd_f32_aaa(dst, x, y, z, length);
  }

  AM_INLINE void abs(float32 *dst, const float32 *x, size_t length) {
    p_abs_f32(dst, x, length);
  }

  AM_INLINE void absCplx(float32 *dst, const float32 *xReal, const float32 *xImag, size_t length) {
    p_absCplx_f32(dst, xReal, xImag, length);
  }

  AM_INLINE void acos(float32 *dst, const float32 *x, size_t length) {
    p_acos_f32(dst, x, length);
  }

  AM_INLINE void asin(float32 *dst, const float32 *x, size_t length) {
    p_asin_f32(dst, x, length);
  }

  AM_INLINE void atan(float32 *dst, const float32 *x, size_t length) {
    p_atan_f32(dst, x, length);
  }

  AM_INLINE void atan2(float32 *dst, const float32 *y, const float32 *x, size_t length) {
    p_atan2_f32(dst, y, x, length);
  }

  AM_INLINE void ceil(float32 *dst, const float32 *x, size_t length) {
    p_ceil_f32(dst, x, length);
  }

  AM_INLINE void cos(float32 *dst, const float32 *x, size_t length) {
    p_cos_f32(dst, x, length);
  }

  AM_INLINE void exp(float32 *dst, const float32 *x, size_t length) {
    p_exp_f32(dst, x, length);
  }

  AM_INLINE void floor(float32 *dst, const float32 *x, size_t length) {
    p_floor_f32(dst, x, length);
  }

  AM_INLINE void log(float32 *dst, const float32 *x, size_t length) {
    p_log_f32(dst, x, length);
  }

  AM_INLINE float32 max(const float32 *x, size_t length) {
    return p_max_f32(x, length);
  }

  AM_INLINE float32 min(const float32 *x, size_t length) {
    return p_min_f32(x, length);
  }

  AM_INLINE void pow(float32 *dst, const float32 *x, float32 y, size_t length) {
    p_pow_f32_as(dst, x, y, length);
  }

  AM_INLINE void pow(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_pow_f32_aa(dst, x, y, length);
  }

  AM_INLINE void random(float32 *dst, float32 low, float32 high, size_t length) {
    m_random->random(dst, low, high, length);
  }

  AM_INLINE void round(float32 *dst, const float32 *x, size_t length) {
    p_round_f32(dst, x, length);
  }

  AM_INLINE void sin(float32 *dst, const float32 *x, size_t length) {
    p_sin_f32(dst, x, length);
  }

  AM_INLINE void sqrt(float32 *dst, const float32 *x, size_t length) {
    p_sqrt_f32(dst, x, length);
  }

  AM_INLINE void tan(float32 *dst, const float32 *x, size_t length) {
    p_tan_f32(dst, x, length);
  }

  AM_INLINE void clamp(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
    p_clamp_f32(dst, x, xMin, xMax, length);
  }

  AM_INLINE void fract(float32 *dst, const float32 *x, size_t length) {
    p_fract_f32(dst, x, length);
  }

  AM_INLINE void ramp(float32 *dst, float32 first, float32 last, size_t length) {
    p_ramp_f32(dst, first, last, length);
  }

  AM_INLINE void sign(float32 *dst, const float32 *x, size_t length) {
    p_sign_f32(dst, x, length);
  }

  AM_INLINE float32 sum(const float32 *x, size_t length) {
    return p_sum_f32(x, length);
  }

  AM_INLINE void sampleLinear(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleLinear_f32(dst, x, t, length, xLength);
  }

  AM_INLINE void sampleLinearRepeat(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleLinearRepeat_f32(dst, x, t, length, xLength);
  }

  AM_INLINE void sampleCubic(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleCubic_f32(dst, x, t, length, xLength);
  }

  AM_INLINE void sampleCubicRepeat(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleCubicRepeat_f32(dst, x, t, length, xLength);
  }

 private:
  // Dispatch table (dynamically assigned function pointers).
  void (*p_add_f32_sa)(float32 *dst, float32 x, const float32 *y, size_t length);
  void (*p_add_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
  void (*p_sub_f32_sa)(float32 *dst, float32 x, const float32 *y, size_t length);
  void (*p_sub_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
  void (*p_mul_f32_sa)(float32 *dst, float32 x, const float32 *y, size_t length);
  void (*p_mul_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
  void (*p_mulCplx_f32_sa)(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length);
  void (*p_mulCplx_f32_aa)(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length);
  void (*p_div_f32_sa)(float32 *dst, float32 x, const float32 *y, size_t length);
  void (*p_div_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
  void (*p_divCplx_f32_sa)(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length);
  void (*p_divCplx_f32_aa)(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length);
  void (*p_madd_f32_saa)(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length);
  void (*p_madd_f32_aaa)(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length);
  void (*p_abs_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_absCplx_f32)(float32 *dst, const float32 *xReal, const float32 *xImag, size_t length);
  void (*p_acos_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_asin_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_atan_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_atan2_f32)(float32 *dst, const float32 *y, const float32 *x, size_t length);
  void (*p_ceil_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_cos_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_exp_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_floor_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_log_f32)(float32 *dst, const float32 *x, size_t length);
  float32 (*p_max_f32)(const float32 *x, size_t length);
  float32 (*p_min_f32)(const float32 *x, size_t length);
  void (*p_pow_f32_as)(float32 *dst, const float32 *x, float32 y, size_t length);
  void (*p_pow_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
  void (*p_round_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_sin_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_sqrt_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_tan_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_clamp_f32)(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length);
  void (*p_fract_f32)(float32 *dst, const float32 *x, size_t length);
  void (*p_ramp_f32)(float32 *dst, float32 first, float32 last, size_t length);
  void (*p_sign_f32)(float32 *dst, const float32 *x, size_t length);
  float32 (*p_sum_f32)(const float32 *x, size_t length);
  void (*p_sampleLinear_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
  void (*p_sampleLinearRepeat_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
  void (*p_sampleCubic_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
  void (*p_sampleCubicRepeat_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);

  // Random number generator.
  Random *m_random;
};

}  // namespace arraymath

#endif // _ARRAYMATH_ARRAYMATH_H
