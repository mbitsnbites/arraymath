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

#ifndef _ARRAYMATH_ARRAYMATHAVX_H
#define _ARRAYMATH_ARRAYMATHAVX_H

#include "common/Architecture.h"

#if defined(AM_USE_X86) && defined(AM_HAS_AVX)

#include "common/Types.h"

namespace arraymath {

class ArrayMathAVX {
 public:
  static void add_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length);
  static void add_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length);
  static void sub_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length);
  static void sub_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length);
  static void mul_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length);
  static void mul_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length);
  static void div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length);
  static void div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length);
  static void madd_f32_saa(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length);
  static void madd_f32_aaa(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length);
  static void sqrt_f32(float32 *dst, const float32 *x, size_t length);
  static void ceil_f32(float32 *dst, const float32 *x, size_t length);
  static void floor_f32(float32 *dst, const float32 *x, size_t length);
  static void round_f32(float32 *dst, const float32 *x, size_t length);
  static void clamp_f32(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length);
  static void fract_f32(float32 *dst, const float32 *x, size_t length);
  static void sign_f32(float32 *dst, const float32 *x, size_t length);
  static void sampleLinear_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
};

}  // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_AVX

#endif // _ARRAYMATH_ARRAYMATHAVX_H
