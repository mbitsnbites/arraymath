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

//------------------------------------------------------------------------------
// This is an generic implementation of the Filter interface that works for all
// CPU architectures. It can handle filters of any order (both FIR and IIR).
//------------------------------------------------------------------------------

#include "FilterGeneric.h"

#include <algorithm>
#include <cstring>

namespace arraymath {

FilterGeneric::FilterGeneric()
    : m_bSize(0), m_aSize(0), m_b(NULL), m_a(NULL), m_bHist(NULL), m_aHist(NULL) {
}

FilterGeneric::~FilterGeneric() {
  if (m_b != NULL)
    delete[] m_b;
  if (m_a != NULL)
    delete[] m_a;
  if (m_bHist != NULL)
    delete[] m_bHist;
  if (m_aHist != NULL)
    delete[] m_aHist;
}

bool FilterGeneric::init(size_t bSize, size_t aSize) {
  // Set filter size.
  m_bSize = bSize;
  m_aSize = aSize;

  // Allocate memory for the coefficients and the history.
  m_b = new float32[bSize];
  m_bHist = new float32[bSize];
  if (!m_b || !m_bHist)
    return false;
  if (aSize > 0) {
    m_a = new float32[aSize];
    m_aHist = new float32[aSize];
    if (!m_a || !m_aHist)
      return false;
  }
  else {
    m_a = m_aHist = NULL;
  }

  // Initialize the filter state.
  std::memset(m_b, 0, bSize * sizeof(float32));
  m_b[0] = 1.0f;
  std::memset(m_bHist, 0, bSize * sizeof(float32));
  if (aSize > 0) {
    std::memset(m_a, 0, aSize * sizeof(float32));
    std::memset(m_aHist, 0, aSize * sizeof(float32));
  }

  return true;
}

void FilterGeneric::setB(const float32 *b) {
  std::memcpy(m_b, b, sizeof(float32) * m_bSize);
}

void FilterGeneric::setA(const float32 *a) {
  std::memcpy(m_a, a, sizeof(float32) * m_aSize);
}

void FilterGeneric::filter(float32 *dst, const float32 *x, size_t length) {
  // Perform run-in part using the history (slow).
  size_t bHistRunIn = m_bSize - 1;
  size_t aHistRunIn = m_aSize;
  size_t k;
  for (k = 0; (bHistRunIn || aHistRunIn) && k < length; ++k) {
    size_t noHistLen, m;

    // FIR part.
    noHistLen = m_bSize - bHistRunIn;
    if (bHistRunIn)
      bHistRunIn--;
    float32 res = m_b[0] * x[k];
    for (m = 1; m < noHistLen; ++m)
      res += m_b[m] * x[k - m];
    for (; m < m_bSize; ++m)
      res += m_b[m] * m_bHist[m - noHistLen];

    // Recursive part.
    noHistLen = m_aSize - aHistRunIn;
    if (aHistRunIn)
      aHistRunIn--;
    for (m = 0; m < noHistLen; ++m)
      res -= m_a[m] * dst[k - 1 - m];
    for (; m < m_aSize; ++m)
      res -= m_a[m] * m_aHist[m - noHistLen];

    dst[k] = res;
  }

  // Perform history-free part (fast).
  for (; k < length; ++k) {
    size_t m;

    // FIR part.
    float32 res = m_b[0] * x[k];
    for (m = 1; m < m_bSize; ++m)
      res += m_b[m] * x[k - m];

    // Recursive part.
    for (m = 0; m < m_aSize; ++m)
      res -= m_a[m] * dst[k - 1 - m];

    dst[k] = res;
  }

  // Update history state.
  size_t histCopy = std::min(m_bSize - 1, length);
  for (k = m_bSize - 2; k >= histCopy; --k)
    m_bHist[k] = m_bHist[k - histCopy];
  for (k = 0; k < histCopy; ++k)
    m_bHist[k] = x[length - 1 - k];
  if (m_aSize > 0) {
    histCopy = std::min(m_aSize, length);
    for (k = m_aSize - 1; k >= histCopy; --k)
      m_aHist[k] = m_aHist[k - histCopy];
    for (k = 0; k < histCopy; ++k)
      m_aHist[k] = dst[length - 1 - k];
  }
}

void FilterGeneric::clearHistory() {
  std::memset(m_bHist, 0, m_bSize * sizeof(float32));
  if (m_aSize > 0) {
    std::memset(m_aHist, 0, m_aSize * sizeof(float32));
  }
}

}  // namespace arraymath
