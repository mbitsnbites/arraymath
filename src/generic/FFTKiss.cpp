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

//------------------------------------------------------------------------------
// This is an generic implementation of the Filter interface that works for all
// CPU architectures. It can handle filters of any order (both FIR and IIR).
//------------------------------------------------------------------------------

#include "generic/FFTKiss.h"

#include <algorithm>
#include <cstring>

namespace arraymath {

FFTKiss::FFTKiss()
    : m_size(0) {
}

FFTKiss::~FFTKiss() {
}

bool FFTKiss::init(size_t size) {
  // TODO(m): Implement me!
  (void)size;
  return true;
}

void FFTKiss::forward(float32* dstReal, float32* dstImag, const float32* x) {
  // TODO(m): Implement me!
  (void)dstReal;
  (void)dstImag;
  (void)x;
}

void FFTKiss::forwardCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // TODO(m): Implement me!
  (void)dstReal;
  (void)dstImag;
  (void)xReal;
  (void)xImag;
}

void FFTKiss::inverse(float32* dst, const float32* xReal, const float32* xImag) {
  // TODO(m): Implement me!
  (void)dst;
  (void)xReal;
  (void)xImag;
}

void FFTKiss::inverseCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // TODO(m): Implement me!
  (void)dstReal;
  (void)dstImag;
  (void)xReal;
  (void)xImag;
}

}  // namespace arraymath
