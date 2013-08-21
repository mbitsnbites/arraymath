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

#ifndef _ARRAYMATH_FILTERFACTORY_H
#define _ARRAYMATH_FILTERFACTORY_H

#include "CPUFeatureDetector.h"
#include "Filter.h"
#include "Types.h"

namespace arraymath {

/// \brief A class for creating filters.
///
/// ### Using the FilterFactory interface
///
/// \code
///    #include <ArrayMath.h>
///    #include <FilterFactory.h>
///
///    void myFunction {
///      // Initialize the FilterFactory object.
///      arraymath::FilterFactory filterFactory;
///
///      // Create a filter object (a 5th order FIR filter).
///      arraymath::Filter* f = filterFactory.createFilter(6, 0);
///      if (!f) {
///        std::cout << "Unable to create filter." << std::endl;
///        return;
///      }
///
///      // Set the filter coefficients.
///      float b[6] = { 0.05f, 0.15f, 0.3f, 0.3f, 0.15f, 0.05f };
///      f->setB(b);
///
///      // Create an input array with white noise.
///      ArrayMath math;
///      const unsigned len = 128;
///      float x[len];
///      math.random(x, -1.0f, 1.0f, len);   // x = random, [-1.0, 1.0)
///
///      // Create a filtered output array.
///      float y[len];
///      f->filter(y, x, len);
///
///      // Delete the filter object.
///      delete f;
///    }
/// \endcode
class FilterFactory {
 public:
  /// Create a new filter object.
  /// \param[in] bSize The number of \c b coefficients.
  /// \param[in] aSize The number of \c a coefficients.
  /// \return A new Filter object, or NULL if none could be created.
  Filter* createFilter(size_t bSize, size_t aSize);

 private:
  CPUFeatureDetector m_cpu;
};

}  // namespace arraymath

#endif // _ARRAYMATH_FILTERFACTORY_H
