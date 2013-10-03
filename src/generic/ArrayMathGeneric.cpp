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

#include "generic/ArrayMathGeneric.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "common/Types.h"

namespace arraymath {

void ArrayMathGeneric::add_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = x + *y++;
  }
}

void ArrayMathGeneric::add_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = *x++ + *y++;
  }
}

void ArrayMathGeneric::sub_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = x - *y++;
  }
}

void ArrayMathGeneric::sub_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = *x++ - *y++;
  }
}

void ArrayMathGeneric::mul_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = x * *y++;
  }
}

void ArrayMathGeneric::mul_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = *x++ * *y++;
  }
}

void ArrayMathGeneric::mulCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  float32 xr = xReal, xi = xImag;
  while (length--) {
    float32 yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathGeneric::mulCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    *dstReal++ = xr * yr - xi * yi;
    *dstImag++ = xr * yi + xi * yr;
  }
}

void ArrayMathGeneric::div_f32_sa(float32 *dst, float32 x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = x / *y++;
  }
}

void ArrayMathGeneric::div_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = *x++ / *y++;
  }
}

void ArrayMathGeneric::divCplx_f32_sa(float32 *dstReal, float32 *dstImag, float32 xReal, float32 xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  float32 xr = xReal, xi = xImag;
  while (length--) {
    float32 yr = *yReal++, yi = *yImag++;
    float32 denom = 1.0f / (yr * yr + yi * yi);
    *dstReal++ = (xr * yr + xi * yi) * denom;
    *dstImag++ = (xi * yr - xr * yi) * denom;
  }
}

void ArrayMathGeneric::divCplx_f32_aa(float32 *dstReal, float32 *dstImag, const float32 *xReal, const float32 *xImag, const float32 *yReal, const float32 *yImag, size_t length) {
  while (length--) {
    float32 xr = *xReal++, xi = *xImag++, yr = *yReal++, yi = *yImag++;
    float32 denom = 1.0f / (yr * yr + yi * yi);
    *dstReal++ = (xr * yr + xi * yi) * denom;
    *dstImag++ = (xi * yr - xr * yi) * denom;
  }
}

void ArrayMathGeneric::madd_f32_saa(float32 *dst, float32 x, const float32 *y, const float32 *z, size_t length) {
  while (length--) {
    *dst++ = x * *y++ + *z++;
  }
}

void ArrayMathGeneric::madd_f32_aaa(float32 *dst, const float32 *x, const float32 *y, const float32 *z, size_t length) {
  while (length--) {
    *dst++ = *x++ + *y++ * *z++;
  }
}

void ArrayMathGeneric::abs_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::fabs(*x++);
  }
}

void ArrayMathGeneric::absCplx_f32(float32 *dst, const float32 *xReal, const float32 *xImag, size_t length) {
  while (length--) {
    float32 re = *xReal++;
    float32 im = *xImag++;
    *dst++ = std::sqrt(re * re + im * im);
  }
}

void ArrayMathGeneric::acos_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::acos(*x++);
  }
}

void ArrayMathGeneric::asin_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::asin(*x++);
  }
}

void ArrayMathGeneric::atan_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::atan(*x++);
  }
}

void ArrayMathGeneric::atan2_f32(float32 *dst, const float32 *y, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::atan2(*y++, *x++);
  }
}

void ArrayMathGeneric::ceil_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::ceil(*x++);
  }
}

void ArrayMathGeneric::cos_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::cos(*x++);
  }
}

void ArrayMathGeneric::exp_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::exp(*x++);
  }
}

void ArrayMathGeneric::floor_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::floor(*x++);
  }
}

void ArrayMathGeneric::log_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::log(*x++);
  }
}

float32 ArrayMathGeneric::max_f32(const float32 *x, size_t length) {
  float32 result = -std::numeric_limits<float>::infinity();
  while (length--) {
    result = std::max(result, *x++);
  }
  return result;
}

float32 ArrayMathGeneric::min_f32(const float32 *x, size_t length) {
  float32 result = std::numeric_limits<float>::infinity();
  while (length--) {
    result = std::min(result, *x++);
  }
  return result;
}

void ArrayMathGeneric::pow_f32_as(float32 *dst, const float32 *x, float32 y, size_t length) {
  while (length--) {
    *dst++ = std::pow(*x++, y);
  }
}

void ArrayMathGeneric::pow_f32_aa(float32 *dst, const float32 *x, const float32 *y, size_t length) {
  while (length--) {
    *dst++ = std::pow(*x++, *y++);
  }
}

void ArrayMathGeneric::round_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::floor(*x++ + 0.5f);
  }
}

void ArrayMathGeneric::sin_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::sin(*x++);
  }
}

void ArrayMathGeneric::sqrt_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::sqrt(*x++);
  }
}

void ArrayMathGeneric::tan_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = std::tan(*x++);
  }
}

void ArrayMathGeneric::clamp_f32(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length) {
  while (length--) {
    float32 val = *x++;
    *dst++ = val < xMin ? xMin : val > xMax ? xMax : val;
  }
}

void ArrayMathGeneric::fract_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    float32 val = *x++;
    *dst++ = val - std::floor(val);
  }
}

void ArrayMathGeneric::fill_f32(float32 *dst, float32 value, size_t length) {
  while (length--) {
    *dst++ = value;
  }
}

void ArrayMathGeneric::ramp_f32(float32 *dst, float32 first, float32 last, size_t length) {
  if (length == 0) {
    return;
  }
  if (length == 1) {
    *dst = first;
    return;
  }

  float32 step = (last - first) / static_cast<float32>(length - 1);
  float32 k = 0.0f;
  while (length--) {
    *dst++ = first + step * k;
    k += 1.0f;
  }
}

void ArrayMathGeneric::sign_f32(float32 *dst, const float32 *x, size_t length) {
  while (length--) {
    *dst++ = *x++ < 0.0f ? -1.0f : 1.0f;
  }
}

float32 ArrayMathGeneric::sum_f32(const float32 *x, size_t length) {
  // TODO(m): We should use pairwise summation or similar here.
  float32 result = 0.0f;
  while (length--) {
    result += *x++;
  }
  return result;
}

void ArrayMathGeneric::sampleLinear_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
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

void ArrayMathGeneric::sampleLinearRepeat_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32 xLengthF = static_cast<float32>(xLength);
  float32 xLengthFInv = 1.0f / xLengthF;
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

void ArrayMathGeneric::sampleCubic_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  while (length--) {
    float32 t2 = *t++;
    t2 = t2 < 0 ? 0 : t2 > maxIdx ? maxIdx : t2;
    size_t idx = std::floor(t2);
    float32 w = t2 - static_cast<float32>(idx);
    float32 w2 = w * w;
    float32 w3 = w2 * w;
    float32 h2 = -2*w3 + 3*w2;
    float32 h1 = 1 - h2;
    float32 h4 = w3 - w2;
    float32 h3 = h4 - w2 + w;
    float32 p1 = x[idx > 0 ? idx - 1 :  0];
    float32 p2 = x[idx];
    float32 p3 = x[idx < maxIdx ? idx + 1 : maxIdx];
    float32 p4 = x[idx < maxIdx - 1 ? idx + 2 : maxIdx];
    *dst++ = h1 * p2 + h2 * p3 + 0.5f * (h3 * (p3 - p1) + h4 * (p4 - p2));
  }
}

void ArrayMathGeneric::sampleCubicRepeat_f32(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength) {
  if (xLength == 0) {
    // If we have nothing to sample, act as if we're sampling only zeros.
    fill_f32(dst, 0.0f, length);
    return;
  }

  size_t maxIdx = xLength - 1;
  float32 xLengthF = static_cast<float32>(xLength);
  float32 xLengthFInv = 1.0f / xLengthF;
  while (length--) {
    float32 t2 = *t++;
    t2 = t2 - std::floor(t2 * xLengthFInv) * xLengthF;
    size_t idx = std::floor(t2);
    float32 w = t2 - static_cast<float32>(idx);
    float32 w2 = w * w;
    float32 w3 = w2 * w;
    float32 h2 = -2*w3 + 3*w2;
    float32 h1 = 1 - h2;
    float32 h4 = w3 - w2;
    float32 h3 = h4 - w2 + w;
    float32 p1 = x[idx > 0 ? idx - 1 : maxIdx];
    float32 p2 = x[idx];
    float32 p3 = x[idx < maxIdx ? idx + 1 : 0];
    float32 p4 = x[idx < maxIdx - 1 ? idx + 2 : idx + 2 - static_cast<size_t>(std::floor((idx + 2) * xLengthFInv) * xLengthF)];
    *dst++ = h1 * p2 + h2 * p3 + 0.5f * (h3 * (p3 - p1) + h4 * (p4 - p2));
  }
}

}  // namespace arraymath
