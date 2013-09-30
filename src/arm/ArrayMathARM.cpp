
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

#include "arm/ArrayMathARM.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "common/Types.h"

namespace arraymath {

namespace {

uint32 asUint32(const float32 x) {
  union {
    float32 f;
    uint32 i;
  } u;
  u.f = x;
  return u.i;
}

}

void ArrayMathARM::abs_f32(float32 *dst, const float32 *x, size_t length) {
  static const uint32 kSignMask = 0x7fffffff;
  uint32 *dst_u32 = reinterpret_cast<uint32*>(dst);
  const uint32 *x_u32 = reinterpret_cast<const uint32*>(x);

  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst_u32) & 15) && length--) {
    *dst_u32++ = (*x_u32++) & kSignMask;
  }

  // 2) Main unrolled loop.
  for (; length >= 4; length -= 4) {
#if defined(__GNUC__)
    __asm__ (
      "ldmia\t%[x_u32]!, {r3-r6}"
      "\n\tbic\tr3, r3, #-2147483648"
      "\n\tbic\tr4, r4, #-2147483648"
      "\n\tbic\tr5, r5, #-2147483648"
      "\n\tbic\tr6, r6, #-2147483648"
      "\n\tstmia\t%[dst_u32]!, {r3-r6}"
      : [dst_u32] "+r" (dst_u32),  [x_u32] "+r" (x_u32)
      :
      : "memory", "r3", "r4", "r5", "r6"
    );
#else
    uint32 x1 = *x_u32++;
    uint32 x2 = *x_u32++;
    uint32 x3 = *x_u32++;
    uint32 x4 = *x_u32++;
    uint32 dst1 = x1 & kSignMask;
    uint32 dst2 = x2 & kSignMask;
    uint32 dst3 = x3 & kSignMask;
    uint32 dst4 = x4 & kSignMask;
    *dst_u32++ = dst1;
    *dst_u32++ = dst2;
    *dst_u32++ = dst3;
    *dst_u32++ = dst4;
#endif
  }

  // 3) Tail loop.
  while (length--) {
    *dst_u32++ = (*x_u32++) & kSignMask;
  }
}

void ArrayMathARM::fill_f32(float32 *dst, float32 value, size_t length) {
  uint32 *dst_u32 = reinterpret_cast<uint32*>(dst);
  uint32 value_u32 = asUint32(value);

  // 1) Align dst to a 16-byte boundary.
  while ((reinterpret_cast<size_t>(dst_u32) & 15) && length--) {
    *dst_u32++ = value_u32;
  }

  // 2) Main unrolled loop.
  for (; length >= 8; length -= 8) {
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
    *dst_u32++ = value_u32;
  }

  // 3) Tail loop.
  while (length--) {
    *dst_u32++ = value_u32;
  }
}

}  // namespace arraymath
