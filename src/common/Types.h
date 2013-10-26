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

#ifndef _ARRAYMATH_TYPES_H
#define _ARRAYMATH_TYPES_H

// Function inline macro.
#if defined(_MSC_VER)
# define AM_INLINE __forceinline
#else
# define AM_INLINE inline
#endif

// Branch optimization macros.
#if defined(__GNUC__)
# define AM_LIKELY(expr) __builtin_expect(!!(expr), 1)
# define AM_UNLIKELY(expr) __builtin_expect(!!(expr), 0)
#else
# define AM_LIKELY(expr) (expr)
# define AM_UNLIKELY(expr) (expr)
#endif

// Needed for size_t.
#include <cstdlib>

namespace arraymath {

typedef unsigned uint32;

typedef float float32;
typedef double float64;

}  // namespace arraymath

#endif // _ARRAYMATH_TYPES_H

