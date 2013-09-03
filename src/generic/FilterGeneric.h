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

#ifndef _ARRAYMATH_FILTERGENERIC_H
#define _ARRAYMATH_FILTERGENERIC_H

#include "Filter.h"

namespace arraymath {

class FilterGeneric : public Filter {
 public:
  virtual ~FilterGeneric();

  virtual void setB(const float32 *b, size_t length);

  virtual void setA(const float32 *a, size_t length);

  virtual void clearHistory();

  virtual void filter(float32 *dst, const float32 *x, size_t length);

 protected:
  FilterGeneric();
  virtual bool init(int bSize, int aSize);

  int runIn(float32 *dst, const float32 *x, size_t length);
  void updateHistory(float32 *dst, const float32 *x, size_t length);

  friend class FilterFactory;

  float32 *m_b;
  float32 *m_a;
  float32 *m_bHist;
  float32 *m_aHist;
  int m_bSize;
  int m_aSize;
};

class FilterGeneric_1_1 : public FilterGeneric {
 public:
  virtual void filter(float32 *dst, const float32 *x, size_t length);

 protected:
  FilterGeneric_1_1() : FilterGeneric() {}

  friend class FilterFactory;
};

class FilterGeneric_3_2 : public FilterGeneric {
 public:
  virtual void filter(float32 *dst, const float32 *x, size_t length);

 protected:
  FilterGeneric_3_2() : FilterGeneric() {}

  friend class FilterFactory;
};

}  // namespace arraymath

#endif // _ARRAYMATH_FILTERGENERIC_H
