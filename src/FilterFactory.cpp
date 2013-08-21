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

#include "FilterFactory.h"
#include "FilterGeneric.h"

namespace arraymath {

Filter* FilterFactory::createFilter(size_t bSize, size_t aSize) {
  // Sanity check.
  if (bSize < 1) {
    return NULL;
  }

  Filter* filter = NULL;
  if (bSize == 3 && aSize == 2) {
    // Optimized filter type: biquad filter.
    filter = new FilterGeneric_3_2();
  }
  else {
    // Generic (fallback) filter that can handle all filter orders.
    filter = new FilterGeneric();
  }

  // Initialize the filter object.
  if (filter) {
    if (!filter->init(bSize, aSize)) {
      delete filter;
      filter = NULL;
    }
  }

  return filter;
}

}  // namespace arraymath
