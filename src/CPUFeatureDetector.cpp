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

#include "CPUFeatureDetector.h"

#ifdef AM_USE_X86

namespace {

enum CPUIDFunc {
  CPUID_VENDOR_ID = 0x00000000,
  CPUID_FEATURES  = 0x00000001
};

void CPUID(CPUIDFunc func, unsigned &a, unsigned &b, unsigned &c, unsigned &d) {
#if defined(__GNUC__)
  __asm__ __volatile__ (
      "cpuid"
      : "=a" (a), "=b" (b), "=c" (c), "=d" (d)
      : "a" (func)
  );
#elif defined(_MSC_VER)
  int info[4];
  __cpuid(info, static_cast<int>(func));
  a = static_cast<unsigned>(info[0]);
  b = static_cast<unsigned>(info[1]);
  c = static_cast<unsigned>(info[2]);
  d = static_cast<unsigned>(info[3]);
#else
  (void)func;
  a = b = c = d = 0;
#endif
}

}  // anonymous namespace

#endif // AM_USE_X86


namespace arraymath {

CPUFeatureDetector::CPUFeatureDetector() {
#ifdef AM_USE_X86
  // Get feature flags.
  unsigned a, b, c, d;
  CPUID(CPUID_VENDOR_ID, a, b, c, d);
  if (a >= CPUID_FEATURES) {
    CPUID(CPUID_FEATURES, a, b, c, d);
  }
  else {
    a = b = c = d = 0;
  }

  // Decode feature flags.
  m_hasMMX = (d & (1 << 23)) != 0;
  m_hasSSE = (d & (1 << 25)) != 0;
  m_hasSSE2 = (d & (1 << 26)) != 0;
  m_hasSSE3 = (c & (1 << 0)) != 0;
  m_hasSSSE3 = (c & (1 << 9)) != 0;
  m_hasSSE4_1 = (c & (1 << 19)) != 0;
  m_hasSSE4_2 = (c & (1 << 20)) != 0;
  m_hasAVX = (c & (1 << 28)) != 0;
  m_hasRDRND = (c & (1 << 30)) != 0;
#endif // AM_USE_X86

#ifdef AM_USE_ARM
  // TODO(m): Check CPU features.
  m_hasNEON = false;
#endif // AM_USE_ARM
}

}  // namespace arraymath
