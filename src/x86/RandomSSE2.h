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

#ifndef _ARRAYMATH_RANDOMSSE2_H
#define _ARRAYMATH_RANDOMSSE2_H

#include "common/Architecture.h"

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)

namespace arraymath {

class Random;

class RandomSSE2Factory {
 public:
  static Random* create();
};

} // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE2

#endif // _ARRAYMATH_RANDOMSSE2_H
