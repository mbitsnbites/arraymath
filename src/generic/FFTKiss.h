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

#ifndef _ARRAYMATH_FFTKISS_H
#define _ARRAYMATH_FFTKISS_H

#ifdef AM_USE_KISS_FFT

#include "FFT.h"

#include "kiss_fft.h"
#include "tools/kiss_fftr.h"

namespace arraymath {

class FFTKiss : public FFT {
 public:
  FFTKiss();
  virtual ~FFTKiss();

  virtual void forward(float32* dstReal, float32* dstImag, const float32* x) const;

  virtual void forwardCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) const;

  virtual void inverse(float32* dst, const float32* xReal, const float32* xImag) const;

  virtual void inverseCplx(float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag) const;

 protected:
  virtual bool init(size_t size);

  friend class FFTFactory;

 private:
  size_t m_size;
  kiss_fft_cfg m_cfg;
  kiss_fft_cfg m_cfgInv;
  kiss_fftr_cfg m_cfgReal;
  kiss_fftr_cfg m_cfgRealInv;
  kiss_fft_cpx* m_inBuf;
  kiss_fft_cpx* m_outBuf;
};

}  // namespace arraymath

#endif // AM_USE_KISS_FFT

#endif // _ARRAYMATH_FFTKISS_H
