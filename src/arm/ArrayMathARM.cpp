
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

#if defined(AM_USE_ARM)

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

  if (length >= 16) {
      // 1) Align x to a 64-byte boundary.
      size_t num_unaligned = reinterpret_cast<size_t>(x_u32) & 63;
      num_unaligned = num_unaligned ? 16 - (num_unaligned >> 2) : 0;
      length -= num_unaligned;
      while (num_unaligned--) {
        *dst_u32++ = (*x_u32++) & kSignMask;
      }

      // 2) Main unrolled loop.
#if defined(__GNUC__)
      __asm__ (
        "\n.abs_loop%=:\n\t"
#if AM_ARM_ARCH >= 6
        "pld   [%[x_u32], #16*4]\n\t"
#endif
        "sub   %[length], %[length], #16\n\t"
        "ldmia %[x_u32]!, {r3,r4,r5,r6,r8,r10,r11,r12}\n\t"
        "bic   r3, r3, #-2147483648\n\t"
        "bic   r4, r4, #-2147483648\n\t"
        "bic   r5, r5, #-2147483648\n\t"
        "bic   r6, r6, #-2147483648\n\t"
        "bic   r8, r8, #-2147483648\n\t"
        "bic   r10, r10, #-2147483648\n\t"
        "bic   r11, r11, #-2147483648\n\t"
        "bic   r12, r12, #-2147483648\n\t"
        "stmia %[dst_u32]!, {r3,r4,r5,r6,r8,r10,r11,r12}\n\t"
        "cmp   %[length], #15\n\t"
        "ldmia %[x_u32]!, {r3,r4,r5,r6,r8,r10,r11,r12}\n\t"
        "bic   r3, r3, #-2147483648\n\t"
        "bic   r4, r4, #-2147483648\n\t"
        "bic   r5, r5, #-2147483648\n\t"
        "bic   r6, r6, #-2147483648\n\t"
        "bic   r8, r8, #-2147483648\n\t"
        "bic   r10, r10, #-2147483648\n\t"
        "bic   r11, r11, #-2147483648\n\t"
        "bic   r12, r12, #-2147483648\n\t"
        "stmia %[dst_u32]!, {r3,r4,r5,r6,r8,r10,r11,r12}\n\t"
        "bhi   .abs_loop%=\n\t"
        : [dst_u32] "+r" (dst_u32), [x_u32] "+r" (x_u32)
        : [length] "r" (length)
        : "memory", "r3", "r4", "r5", "r6", "r8", "r10", "r11", "r12"
      );
#endif
      for (; length >= 8; length -= 8) {
        uint32 x1 = *x_u32++;
        uint32 x2 = *x_u32++;
        uint32 x3 = *x_u32++;
        uint32 x4 = *x_u32++;
        uint32 x5 = *x_u32++;
        uint32 x6 = *x_u32++;
        uint32 x7 = *x_u32++;
        uint32 x8 = *x_u32++;
        *dst_u32++ = x1 & kSignMask;
        *dst_u32++ = x2 & kSignMask;
        *dst_u32++ = x3 & kSignMask;
        *dst_u32++ = x4 & kSignMask;
        *dst_u32++ = x5 & kSignMask;
        *dst_u32++ = x6 & kSignMask;
        *dst_u32++ = x7 & kSignMask;
        *dst_u32++ = x8 & kSignMask;
      }
  }

  // 3) Tail loop.
  while (length--) {
    *dst_u32++ = (*x_u32++) & kSignMask;
  }
}

void ArrayMathARM::fill_f32(float32 *dst, float32 value, size_t length) {
  uint32 *dst_u32 = reinterpret_cast<uint32*>(dst);
  uint32 value_u32 = asUint32(value);

  if (length >= 64) {
      // 1) Align dst to a 64-byte boundary.
      size_t num_unaligned = reinterpret_cast<size_t>(dst_u32) & 63;
      num_unaligned = num_unaligned ? 16 - (num_unaligned >> 2) : 0;
      length -= num_unaligned;
      while (num_unaligned--) {
        *dst_u32++ = value_u32;
      }

      // 2) Main unrolled loop.
#if defined(__GNUC__)
      __asm__ (
        "mov   r3, %[value_u32]\n\t"
        "mov   r4, %[value_u32]\n\t"
        "mov   r5, %[value_u32]\n\t"
        "mov   r6, %[value_u32]\n\t"
        "mov   r8, %[value_u32]\n\t"
        "mov   r10, %[value_u32]\n\t"
        "mov   r11, %[value_u32]\n\t"
        "mov   r12, %[value_u32]\n"
        ".fill_loop%=:\n\t"
        "sub   %[length], %[length], #8\n\t"
        "stmia %[dst_u32]!, {r3,r4,r5,r6,r8,r10,r11,r12}\n\t"
        "cmp   %[length], #7\n\t"
        "bhi   .fill_loop%=\n\t"
        : [dst_u32] "+r" (dst_u32)
        : [value_u32] "r" (value_u32), [length] "r" (length)
        : "memory", "r3", "r4", "r5", "r6", "r8", "r10", "r11", "r12"
      );
#else
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
#endif
  }

  // 3) Tail loop.
  while (length--) {
    *dst_u32++ = value_u32;
  }
}

}  // namespace arraymath

#endif // AM_USE_ARM
