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

#include "Types.h"

namespace arraymath {

class CPUFeatureDetector {
 public:
  CPUFeatureDetector();

#ifdef AM_USE_X86
  bool hasMMX() const {
    return m_hasMMX;
  }

  bool hasSSE() const {
    return m_hasSSE;
  }

  bool hasSSE2() const {
    return m_hasSSE2;
  }

  bool hasSSE3() const {
    return m_hasSSE3;
  }

  bool hasSSSE3() const {
    return m_hasSSSE3;
  }

  bool hasSSE4_1() const {
    return m_hasSSE4_1;
  }

  bool hasSSE4_2() const {
    return m_hasSSE4_2;
  }

  bool hasAVX() const {
    return m_hasAVX;
  }

  bool hasRDRND() const {
    return m_hasRDRND;
  }
#endif // AM_USE_X86

#ifdef AM_USE_ARM
  bool hasNEON() const {
    return m_hasNEON;
  }
#endif // AM_USE_ARM

 private:
#ifdef AM_USE_X86
  bool m_hasMMX;
  bool m_hasSSE;
  bool m_hasSSE2;
  bool m_hasSSE3;
  bool m_hasSSSE3;
  bool m_hasSSE4_1;
  bool m_hasSSE4_2;
  bool m_hasAVX;
  bool m_hasRDRND;
#endif // AM_USE_X86

#ifdef AM_USE_ARM
  bool m_hasNEON;
#endif // AM_USE_ARM
};

}  // namespace arraymath

#endif // _ARRAYMATH_CPUFEATUREDETECTOR_H
