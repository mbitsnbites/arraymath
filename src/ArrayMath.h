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
#include "common/Types.h"

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

  /// Complex multiplication.
  /// Compute (dstReal + j*dstImag) = (xReal + j*xImag) * (yReal + j*yImag).
  /// \param[out] dstReal Destination array, real part (\c length elements).
  /// \param[out] dstImag Destination array, imaginary part (\c length elements).
  /// \param[in] xReal Source scalar, real part.
  /// \param[in] xImag Source scalar, imaginary part.
  /// \param[in] yReal Source array, real part (\c length elements).
  /// \param[in] yImag Source array, imaginary part (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void mulCplx(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_mulCplx_f32_sa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  /// Complex multiplication.
  /// Compute (dstReal + j*dstImag) = (xReal + j*xImag) * (yReal + j*yImag).
  /// \param[out] dstReal Destination array, real part (\c length elements).
  /// \param[out] dstImag Destination array, imaginary part (\c length elements).
  /// \param[in] xReal Source array, real part (\c length elements).
  /// \param[in] xImag Source array, imaginary part (\c length elements).
  /// \param[in] yReal Source array, real part (\c length elements).
  /// \param[in] yImag Source array, imaginary part (\c length elements).
  /// \param[in] length Number of elements to process.
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

  /// Complex division.
  /// Compute (dstReal + j*dstImag) = (xReal + j*xImag) / (yReal + j*yImag).
  /// \param[out] dstReal Destination array, real part (\c length elements).
  /// \param[out] dstImag Destination array, imaginary part (\c length elements).
  /// \param[in] xReal Source scalar, real part.
  /// \param[in] xImag Source scalar, imaginary part.
  /// \param[in] yReal Source array, real part (\c length elements).
  /// \param[in] yImag Source array, imaginary part (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void divCplx(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
    p_divCplx_f32_sa(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }

  /// Complex division.
  /// Compute (dstReal + j*dstImag) = (xReal + j*xImag) / (yReal + j*yImag).
  /// \param[out] dstReal Destination array, real part (\c length elements).
  /// \param[out] dstImag Destination array, imaginary part (\c length elements).
  /// \param[in] xReal Source array, real part (\c length elements).
  /// \param[in] xImag Source array, imaginary part (\c length elements).
  /// \param[in] yReal Source array, real part (\c length elements).
  /// \param[in] yImag Source array, imaginary part (\c length elements).
  /// \param[in] length Number of elements to process.
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

  /// Compute abs(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void abs(float32 *dst, const float32 *x, size_t length) {
    p_abs_f32(dst, x, length);
  }

  /// Complex absolute value.
  /// Compute dst = |xReal + j*xImag|
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] xReal Source array, real part (\c length elements).
  /// \param[in] xImag Source array, imaginary part (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void absCplx(float32 *dst, const float32 *xReal, const float32 *xImag, size_t length) {
    p_absCplx_f32(dst, xReal, xImag, length);
  }

  /// Compute acos(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void acos(float32 *dst, const float32 *x, size_t length) {
    p_acos_f32(dst, x, length);
  }

  /// Compute asin(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void asin(float32 *dst, const float32 *x, size_t length) {
    p_asin_f32(dst, x, length);
  }

  /// Compute atan(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void atan(float32 *dst, const float32 *x, size_t length) {
    p_atan_f32(dst, x, length);
  }

  /// Compute atan2(y, x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void atan2(float32 *dst, const float32 *y, const float32 *x, size_t length) {
    p_atan2_f32(dst, y, x, length);
  }

  /// Compute ceil(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void ceil(float32 *dst, const float32 *x, size_t length) {
    p_ceil_f32(dst, x, length);
  }

  /// Compute cos(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void cos(float32 *dst, const float32 *x, size_t length) {
    p_cos_f32(dst, x, length);
  }

  /// Compute exp(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void exp(float32 *dst, const float32 *x, size_t length) {
    p_exp_f32(dst, x, length);
  }

  /// Compute floor(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void floor(float32 *dst, const float32 *x, size_t length) {
    p_floor_f32(dst, x, length);
  }

  /// Compute log(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void log(float32 *dst, const float32 *x, size_t length) {
    p_log_f32(dst, x, length);
  }

  /// Find the maximum value in the source array.
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \return The maxium value found in the array \c x.
  AM_INLINE float32 max(const float32 *x, size_t length) {
    return p_max_f32(x, length);
  }

  /// Find the minimum value in the source array.
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \return The minimum value found in the array \c x.
  AM_INLINE float32 min(const float32 *x, size_t length) {
    return p_min_f32(x, length);
  }

  /// Compute pow(x, y) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source scalar.
  /// \param[in] length Number of elements to process.
  AM_INLINE void pow(float32 *dst, const float32 *x, float32 y, size_t length) {
    p_pow_f32_as(dst, x, y, length);
  }

  /// Compute pow(x, y) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] y Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void pow(float32 *dst, const float32 *x, const float32 *y, size_t length) {
    p_pow_f32_aa(dst, x, y, length);
  }

  /// Generate random numbers and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] low Lowest value to produce (inclusive).
  /// \param[in] high Highest value to produce (exclusive).
  /// \param[in] length Number of elements to process.
  AM_INLINE void random(float32 *dst, float32 low, float32 high, size_t length) {
    m_random_f32->random(dst, low, high, length);
  }

  /// Compute round(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void round(float32 *dst, const float32 *x, size_t length) {
    p_round_f32(dst, x, length);
  }

  /// Compute sin(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void sin(float32 *dst, const float32 *x, size_t length) {
    p_sin_f32(dst, x, length);
  }

  /// Compute sqrt(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void sqrt(float32 *dst, const float32 *x, size_t length) {
    p_sqrt_f32(dst, x, length);
  }

  /// Compute tan(x) and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void tan(float32 *dst, const float32 *x, size_t length) {
    p_tan_f32(dst, x, length);
  }

  /// Clamp the input array and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] xMin Minimum value to clamp to.
  /// \param[in] xMax Maximum value to clamp to.
  /// \param[in] length Number of elements to process.
  AM_INLINE void clamp(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
    p_clamp_f32(dst, x, xMin, xMax, length);
  }

  /// Extract the fractional part of the input array and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void fract(float32 *dst, const float32 *x, size_t length) {
    p_fract_f32(dst, x, length);
  }

  /// Generate a linear ramp.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] first First value.
  /// \param[in] last Last value.
  /// \param[in] length Number of elements to process.
  AM_INLINE void ramp(float32 *dst, float32 first, float32 last, size_t length) {
    p_ramp_f32(dst, first, last, length);
  }

  /// Extract the sign of the input array and store the result in dst.
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  AM_INLINE void sign(float32 *dst, const float32 *x, size_t length) {
    p_sign_f32(dst, x, length);
  }

  /// Calculate the sum of all elements in an array.
  /// \param[in] x Source array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \return The sum of all elements in the array \c x.
  AM_INLINE float32 sum(const float32 *x, size_t length) {
    return p_sum_f32(x, length);
  }

  /// Perform linear interpolation (clamping index).
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Value array (\c xLength elements).
  /// \param[in] t Sample positions array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \param[in] xLength Number of elements in the value array.
  AM_INLINE void sampleLinear(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleLinear_f32(dst, x, t, length, xLength);
  }

  /// Perform linear interpolation (repeating index).
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Value array (\c xLength elements).
  /// \param[in] t Sample positions array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \param[in] xLength Number of elements in the value array.
  AM_INLINE void sampleLinearRepeat(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleLinearRepeat_f32(dst, x, t, length, xLength);
  }

  /// Perform cubic interpolation (clamping index).
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Value array (\c xLength elements).
  /// \param[in] t Sample positions array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \param[in] xLength Number of elements in the value array.
  AM_INLINE void sampleCubic(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
    p_sampleCubic_f32(dst, x, t, length, xLength);
  }

  /// Perform cubic interpolation (repeating index).
  /// \param[out] dst Destination array (\c length elements).
  /// \param[in] x Value array (\c xLength elements).
  /// \param[in] t Sample positions array (\c length elements).
  /// \param[in] length Number of elements to process.
  /// \param[in] xLength Number of elements in the value array.
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
  Random *m_random_f32;
};

}  // namespace arraymath

#endif // _ARRAYMATH_ARRAYMATH_H
