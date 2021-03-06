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

#ifndef _ARRAYMATH_CPUFEATUREDETECTOR_H
#define _ARRAYMATH_CPUFEATUREDETECTOR_H

#include "common/Architecture.h"

namespace arraymath {

class CPUFeatureDetector {
 public:
  CPUFeatureDetector();

#ifdef AM_USE_X86
  bool hasSSE() const {
    return m_hasSSE;
  }

  bool hasSSE2() const {
    return m_hasSSE2;
  }

  bool hasSSE4() const {
    return m_hasSSE4;
  }

  bool hasAVX() const {
    return m_hasAVX;
  }
#endif // AM_USE_X86

#ifdef AM_USE_ARM
  bool hasNEON() const {
    return m_hasNEON;
  }

  bool hasNEON_FMA() const {
    return m_hasNEON_FMA;
  }
#endif // AM_USE_ARM

 private:
#ifdef AM_USE_X86
  bool m_hasSSE;
  bool m_hasSSE2;
  bool m_hasSSE4;
  bool m_hasAVX;
#endif // AM_USE_X86

#ifdef AM_USE_ARM
  bool m_hasNEON;
  bool m_hasNEON_FMA;
#endif // AM_USE_ARM
};

}  // namespace arraymath

#endif // _ARRAYMATH_CPUFEATUREDETECTOR_H
