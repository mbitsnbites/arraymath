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

#include "x86/FilterAVX.h"

#include <algorithm>
#include <cstring>

#include <immintrin.h>

namespace arraymath {

namespace {

AM_INLINE
float horizontalSum(__m256 x) {
#if 1
  // From http://stackoverflow.com/questions/13219146/how-to-sum-m256-horizontally
  const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
  const __m128 loQuad = _mm256_castps256_ps128(x);
  const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
  const __m128 loDual = sumQuad;
  const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
  const __m128 sumDual = _mm_add_ps(loDual, hiDual);
  const __m128 lo = sumDual;
  const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
  const __m128 sum = _mm_add_ss(lo, hi);
  return _mm_cvtss_f32(sum);
#else
  // This is slightly slower on my Core i7...
  const __m128 xHi = _mm256_extractf128_ps(x, 1);
  const __m128 xLo = _mm256_castps256_ps128(x);
  __m128 sum = _mm_add_ps(xLo, xHi);
  sum = _mm_hadd_ps(sum, sum);
  sum = _mm_hadd_ps(sum, sum);
  return _mm_cvtss_f32(sum);
#endif
}

}

FilterAVX_FIR::FilterAVX_FIR()
    : m_b(NULL), m_history(NULL), m_size(0) {
}

FilterAVX_FIR::~FilterAVX_FIR() {
  if (m_b != NULL)
    delete[] m_b;
  if (m_history != NULL)
    delete[] m_history;
}

bool FilterAVX_FIR::init(int size) {
  // Set filter size.
  m_size = size;

  // Allocate memory for the coefficients.
  // NOTE: Make sure that the buffer is 32-byte aligned!
  int numAVXElements = (size + 7) / 8;
  m_b = reinterpret_cast<float32*>(new __m256[numAVXElements]);
  if (!m_b)
    return false;

  // Allocate memory for the history + additional space for concatenating the
  // first bSize - 1 samples of the input buffer to the history during the
  // run-in step.
  // NOTE: The buffer does not have to be AVX-aligned.
  m_history = new float32[(size - 1) * 2];
  if (!m_history)
    return false;

  // Initialize the filter state.
  std::memset(m_b, 0, (size - 1) * sizeof(float32));
  m_b[size - 1] = 1.0f;
  std::memset(m_history, 0, (size - 1) * sizeof(float32));

  return true;
}

void FilterAVX_FIR::setB(const float32 *b, size_t length) {
  int len = std::min(m_size, static_cast<int>(length));

  // NOTE: The impulse response is reversed in order to enable incremental
  // memory access during the filtering step.
  float32* dst = &m_b[m_size - 1];
  while (len--)
    *dst-- = *b++;
}

void FilterAVX_FIR::clearHistory() {
  std::memset(m_history, 0, (m_size - 1) * sizeof(float32));
}

void FilterAVX_FIR::filter(float32 *dst, const float32 *x, size_t length) {
  const int len = static_cast<int>(length);
  const int histLen = m_size - 1;
  const int runInLen = std::min(len, histLen);

  // Append the required number of samples from the input buffer to the history
  // buffer.
  std::memcpy(&m_history[histLen], x, runInLen * sizeof(float32));

  // Perform run-in part using the history.
  const float32 *in = m_history;
  float32 *out = dst;
  for (int k = 0; k < runInLen; ++k) {
    const float32 *src = in;
    const float32 *b = m_b;
    int bSize = m_size;
    __m256 _y = _mm256_setzero_ps();
    for (; bSize >= 16; bSize -= 16) {
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b), _mm256_loadu_ps(src)));
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b + 8), _mm256_loadu_ps(src + 8)));
      b += 16; src += 16;
    }
    for (; bSize >= 8; bSize -= 8) {
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b), _mm256_loadu_ps(src)));
      b += 8; src += 8;
    }
    float32 y = horizontalSum(_y);
    while (bSize--)
      y += *b++ * *src++;
    *out++ = y;
    ++in;
  }

  int samplesLeft = len - runInLen;

  // Main history-less loop.
  in = &x[runInLen];
  while (samplesLeft--) {
    const float32 *src = in - (m_size - 1);
    const float32 *b = m_b;
    int bSize = m_size;
    __m256 _y = _mm256_setzero_ps();
    for (; bSize >= 16; bSize -= 16) {
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b), _mm256_loadu_ps(src)));
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b + 8), _mm256_loadu_ps(src + 8)));
      b += 16; src += 16;
    }
    for (; bSize >= 8; bSize -= 8) {
      _y = _mm256_add_ps(_y, _mm256_mul_ps(_mm256_load_ps(b), _mm256_loadu_ps(src)));
      b += 8; src += 8;
    }
    float32 y = horizontalSum(_y);
    while (bSize--)
      y += *b++ * *src++;
    *out++ = y;
    ++in;
  }

  // Shift history.
  const int histCopy = std::max((m_size - 1) - len, 0);
  if (histCopy > 0) {
    const float32 *src = &m_history[m_size - 2];
    float32 *dst = &m_history[histCopy - 1];
    std::memmove(dst, src, histCopy * sizeof(float32));
  }

  // Copy from the tail of x to the history.
  const int xCopy = (m_size - 1) - histCopy;
  if (xCopy > 0) {
    const float32 *src = &x[len - xCopy];
    float32 *dst = &m_history[histCopy];
    std::memcpy(dst, src, xCopy * sizeof(float32));
  }
}

}  // namespace arraymath
