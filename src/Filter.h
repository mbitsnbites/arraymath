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

#ifndef _ARRAYMATH_FILTER_H
#define _ARRAYMATH_FILTER_H

#include "common/Types.h"

namespace arraymath {

/// \brief Filter interface.
class Filter {
 public:
  virtual ~Filter() {}

  /// Set the \c b coefficients of the filter.
  /// \param[in] b The \c b coefficients.
  /// \note The number of elements in the array \c b must be equal to the
  /// \c bSize argument of the filter, as given when creating the filter.
  virtual void setB(const float32 *b) = 0;

  /// Set the \c a coefficients of the filter.
  /// \param[in] a The \c a coefficients.
  /// \note The number of elements in the array \c a must be equal to the
  /// \c aSize argument of the filter, as given when creating the filter.
  virtual void setA(const float32 *a) = 0;

  /// Clear the history state of the filter.
  virtual void clearHistory() = 0;

  /// Filter the given array.
  /// \param[out] dst The destination array (\c length elements).
  /// \param[in] x The source array (\c length elements).
  /// \param[in] length Number of elements to process.
  virtual void filter(float32 *dst, const float32 *x, size_t length) = 0;

 protected:
  virtual bool init(int bSize, int aSize) = 0;

  friend class FilterFactory;
};

}  // namespace arraymath

#endif // _ARRAYMATH_FILTER_H
