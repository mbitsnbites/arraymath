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

#ifndef _ARRAYMATH_FILTERAVX_H
#define _ARRAYMATH_FILTERAVX_H

#include "Filter.h"

namespace arraymath {

class FilterAVX_FIR : public Filter {
 public:
  virtual ~FilterAVX_FIR();
  virtual void setB(const float32 *b, size_t length);
  virtual void setA(const float32*, size_t) {}
  virtual void clearHistory();
  virtual void filter(float32 *dst, const float32 *x, size_t length);

 protected:
  FilterAVX_FIR();

  bool init(int size);

  virtual bool init(int bSize, int) {
    return init(bSize);
  }

  friend class FilterFactory;

  float32 *m_b;
  float32 *m_history;
  int m_size;
};

}  // namespace arraymath

#endif // _ARRAYMATH_FILTERAVX_H
