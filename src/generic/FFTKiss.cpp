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
// This is an interface to the Kiss FFT library.
// TODO(m): It's inefficient to have to convert between de-interleaved and
// interleaved data representations twice (!) for each transformation call.
//------------------------------------------------------------------------------

#ifdef AM_USE_KISS_FFT

#include "generic/FFTKiss.h"

namespace arraymath {

FFTKiss::FFTKiss()
    : m_size(0), m_cfg(NULL), m_cfgInv(NULL), m_inBuf(NULL), m_outBuf(NULL) {
}

FFTKiss::~FFTKiss() {
  if (m_cfg) {
    kiss_fft_free(m_cfg);
  }
  if (m_cfgInv) {
    kiss_fft_free(m_cfgInv);
  }
  if (m_inBuf) {
    delete[] m_inBuf;
  }
  if (m_outBuf) {
    delete[] m_outBuf;
  }
}

bool FFTKiss::init(size_t size) {
  m_size = size;

  // Allocate Kiss FFT state.
  m_cfg = kiss_fft_alloc(size, 0, NULL, NULL);
  if (!m_cfg) {
    return false;
  }
  m_cfgInv = kiss_fft_alloc(size, 1, NULL, NULL);
  if (!m_cfgInv) {
    return false;
  }

  // TODO(m): How about real valued configurations?

  // Allocate memory for temporary working buffers.
  m_inBuf = new kiss_fft_cpx[size];
  if (!m_inBuf) {
    return false;
  }
  m_outBuf = new kiss_fft_cpx[size];
  if (!m_outBuf) {
    return false;
  }

  return true;
}

void FFTKiss::forward(float32* dstReal, float32* dstImag, const float32* x) {
  // TODO(m): Implement me!
  (void)dstReal;
  (void)dstImag;
  (void)x;
}

void FFTKiss::forwardCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // Copy input signal to an interleaved buffer.
  for (size_t i = 0; i < m_size; ++i) {
    m_inBuf[i].r = *xReal++;
    m_inBuf[i].i = *xImag++;
  }

  // Perform forward transform.
  kiss_fft(m_cfg, m_inBuf, m_outBuf);

  // Copy interleaved output buffer to the output signal.
  for (size_t i = 0; i < m_size; ++i) {
    *dstReal++ = m_outBuf[i].r;
    *dstImag++ = m_outBuf[i].i;
  }
}

void FFTKiss::inverse(float32* dst, const float32* xReal, const float32* xImag) {
  // TODO(m): Implement me!
  (void)dst;
  (void)xReal;
  (void)xImag;
}

void FFTKiss::inverseCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // Copy input signal to an interleaved buffer.
  for (size_t i = 0; i < m_size; ++i) {
    m_inBuf[i].r = *xReal++;
    m_inBuf[i].i = *xImag++;
  }

  // Perform forward transform.
  kiss_fft(m_cfgInv, m_inBuf, m_outBuf);

  // Copy interleaved output buffer to the output signal.
  for (size_t i = 0; i < m_size; ++i) {
    *dstReal++ = m_outBuf[i].r;
    *dstImag++ = m_outBuf[i].i;
  }
}

}  // namespace arraymath

#endif // AM_USE_KISS_FFT
