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
//
// TODO(m): Check / define the order of the frequency components.
//
// TODO(m): It's inefficient to have to convert between de-interleaved and
// interleaved data representations for each transformation call.
//
// TODO(m): Could we apply the scaling factor directly in the Kiss FFT routines
// instead?
//------------------------------------------------------------------------------

#ifdef AM_USE_KISS_FFT

#include "generic/FFTKiss.h"

#include <cmath>

namespace arraymath {

FFTKiss::FFTKiss()
    : m_size(0), m_cfg(NULL), m_cfgInv(NULL), m_cfgReal(NULL), m_cfgRealInv(NULL), m_inBuf(NULL), m_outBuf(NULL) {
}

FFTKiss::~FFTKiss() {
  if (m_cfg) {
    kiss_fft_free(m_cfg);
  }
  if (m_cfgInv) {
    kiss_fft_free(m_cfgInv);
  }
  if (m_cfgReal) {
    kiss_fftr_free(m_cfgReal);
  }
  if (m_cfgRealInv) {
    kiss_fftr_free(m_cfgRealInv);
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
  if ((size & 1) == 0) {
    // kiss_fftr_alloc only allows even sizes.
    m_cfgReal = kiss_fftr_alloc(size, 0, NULL, NULL);
    if (!m_cfgReal) {
      return false;
    }
    m_cfgRealInv = kiss_fftr_alloc(size, 1, NULL, NULL);
    if (!m_cfgRealInv) {
      return false;
    }
  }

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
  // Scaling factor.
  const float32 scale = 1.0f / std::sqrt(static_cast<float>(m_size));

  // If size is odd, we need to use the complex version.
  if (m_size & 1) {
    // Copy input signal to an interleaved buffer.
    for (size_t i = 0; i < m_size; ++i) {
      m_inBuf[i].r = scale * *x++;
      m_inBuf[i].i = 0.0f;
    }

    // Perform forward transform.
    kiss_fft(m_cfg, m_inBuf, m_outBuf);

    // Copy interleaved output buffer to the output signal.
    for (size_t i = 0; i < m_size; ++i) {
      *dstReal++ = m_outBuf[i].r;
      *dstImag++ = m_outBuf[i].i;
    }
  }
  else {
    // Perform forward transform.
    kiss_fftr(m_cfgReal, x, m_outBuf);

    // Copy interleaved output buffer to the output signal.
    // Note: Here we duplicate / mirror the lower part of the spectrum, in order
    // to mimic the result of doing a regular complex FFT of a signal with the
    // imaginary part set to zero.
    size_t numOutSamples = (m_size / 2) + 1;
    for (size_t i = 0; i < numOutSamples; ++i) {
      dstReal[i] = dstReal[(m_size - 1) - i] = scale * m_outBuf[i].r;
      dstImag[i] = dstImag[(m_size - 1) - i] = scale * m_outBuf[i].i;
    }
  }
}

void FFTKiss::forwardCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // Scaling factor.
  const float32 scale = 1.0f / std::sqrt(static_cast<float>(m_size));

  // Copy input signal to an interleaved buffer.
  for (size_t i = 0; i < m_size; ++i) {
    m_inBuf[i].r = *xReal++;
    m_inBuf[i].i = *xImag++;
  }

  // Perform forward transform.
  kiss_fft(m_cfg, m_inBuf, m_outBuf);

  // Copy interleaved output buffer to the output signal.
  for (size_t i = 0; i < m_size; ++i) {
    *dstReal++ = scale * m_outBuf[i].r;
    *dstImag++ = scale * m_outBuf[i].i;
  }
}

void FFTKiss::inverse(float32* dst, const float32* xReal, const float32* xImag) {
  // Scaling factor.
  const float32 scale = 1.0f / std::sqrt(static_cast<float>(m_size));

  // If size is odd, we need to use the complex version.
  if (m_size & 1) {
    // Copy input signal to an interleaved buffer.
    for (size_t i = 0; i < m_size; ++i) {
      m_inBuf[i].r = *xReal++;
      m_inBuf[i].i = *xImag++;
    }

    // Perform reverse transform.
    kiss_fft(m_cfgInv, m_inBuf, m_outBuf);

    // Copy interleaved output buffer to the output signal.
    for (size_t i = 0; i < m_size; ++i) {
      *dst++ = scale * m_outBuf[i].r;
    }
  }
  else {
    // Copy input signal to an interleaved buffer.
    // Note: We only use the lower half of the spectrum, since kiss_fftri will
    // not use the upper part (it's symmetric).
    size_t numInSamples = (m_size / 2) + 1;
    for (size_t i = 0; i < numInSamples; ++i) {
      m_inBuf[i].r = scale * *xReal++;
      m_inBuf[i].i = scale * *xImag++;
    }

    // Perform inverse transform.
    kiss_fftri(m_cfgRealInv, m_inBuf, dst);
  }
}

void FFTKiss::inverseCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) {
  // Scaling factor.
  const float32 scale = 1.0f / std::sqrt(static_cast<float>(m_size));

  // Copy input signal to an interleaved buffer.
  for (size_t i = 0; i < m_size; ++i) {
    m_inBuf[i].r = *xReal++;
    m_inBuf[i].i = *xImag++;
  }

  // Perform forward transform.
  kiss_fft(m_cfgInv, m_inBuf, m_outBuf);

  // Copy interleaved output buffer to the output signal.
  for (size_t i = 0; i < m_size; ++i) {
    *dstReal++ = scale * m_outBuf[i].r;
    *dstImag++ = scale * m_outBuf[i].i;
  }
}

}  // namespace arraymath

#endif // AM_USE_KISS_FFT
