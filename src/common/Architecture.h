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

#ifndef _ARRAYMATH_ARCHITECTURE_H
#define _ARRAYMATH_ARCHITECTURE_H

// Operating system detection.
#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
# define AM_OS_WINDOWS
#elif defined(__APPLE__)
# include <TargetConditionals.h>
# if defined(TARGET_OS_IPHONE) || defined(TARGET_OS_EMBEDDED) || defined(TARGET_IPHONE_SIMULATOR)
#  define AM_OS_IOS
# elif defined(TARGET_OS_MAC)
#  define AM_OS_MACOSX
# endif
#elif defined(__ANDROID__)
# define AM_OS_ANDROID
#elif defined(unix) || defined(__unix__) || defined(__unix) || defined(__linux__)
# define AM_OS_UNIX
#endif

// CPU architecture detection.
#if defined(_M_IX86) || defined(__i386__) || defined(_X86_) || defined(_M_X64) || defined(__amd64__)

// We're compiling for an x86 machine.
# define AM_USE_X86 1

// x86 compile time feature detection.
// NOTE: You should set these manually in your build config / project.
# if !defined(AM_HAS_SSE) && (defined(__SSE__) || defined(_M_X64))
#  define AM_HAS_SSE 1
# endif
# if !defined(AM_HAS_SSE2) && (defined(__SSE2__) || defined(_M_X64))
#  define AM_HAS_SSE2 1
# endif
# if !defined(AM_HAS_AVX) && defined(__AVX__)
#  define AM_HAS_AVX 1
# endif

#elif defined(_M_ARM) || defined(__arm__)

// We're compiling for an ARM machine.
# define AM_USE_ARM 1

// ARM compile time feature detection.
// NOTE: You should set these manually in your build config / project.
# if !defined(AM_HAS_NEON) && defined(__ARM_NEON__)
#  define AM_HAS_NEON 1
# endif

#endif


#endif // _ARRAYMATH_ARCHITECTURE_H
