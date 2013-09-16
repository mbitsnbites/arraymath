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

#include "common/CPUFeatureDetector.h"

#ifdef AM_USE_X86

#if defined(__GNUC__) || defined(__clang__)
#include <cpuid.h>
#endif
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace {

enum CPUIDFunc {
  CPUID_VENDOR_ID = 0x00000000,
  CPUID_FEATURES  = 0x00000001
};

void CPUID(CPUIDFunc func, unsigned &a, unsigned &b, unsigned &c, unsigned &d) {
#if defined(__GNUC__) || defined(__clang__)
  __get_cpuid(func, &a, &b, &c, &d);
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

#if defined(AM_USE_ARM) && defined(AM_OS_ANDROID)
#include <cpu-features.h>
#endif // AM_USE_ARM && AM_OS_ANDROID


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
//  m_hasMMX = (d & (1 << 23)) != 0;
  m_hasSSE = (d & (1 << 25)) != 0;
  m_hasSSE2 = (d & (1 << 26)) != 0;
//  m_hasSSE3 = (c & (1 << 0)) != 0;
//  m_hasSSSE3 = (c & (1 << 9)) != 0;
//  m_hasSSE4_1 = (c & (1 << 19)) != 0;
//  m_hasSSE4_2 = (c & (1 << 20)) != 0;
  m_hasAVX = (c & (1 << 28)) != 0;
//  m_hasRDRND = (c & (1 << 30)) != 0;
#endif // AM_USE_X86

#ifdef AM_USE_ARM
# if defined(AM_OS_ANDROID)
  // Get feature flags.
  uint64_t features = android_getCpuFeatures();

  // Decode feature flags.
  m_hasNEON = (features & ANDROID_CPU_ARM_FEATURE_NEON) != 0;
  m_hasNEON_FMA = (features & ANDROID_CPU_ARM_FEATURE_NEON_FMA) != 0;
# else
  // Fall back to compile-time detection. Can we do better?
#  ifdef AM_HAS_NEON
  m_hasNEON = true;
  m_hasNEON_FMA = false;  // TODO(m): Can we do better?
#  else
  m_hasNEON = false;
  m_hasNEON_FMA = false;
#  endif // AM_HAS_NEON
# endif // AM_OS_ANDROID
#endif // AM_USE_ARM
}

}  // namespace arraymath
