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

#ifndef _ARRAYMATH_FFTFACTORY_H
#define _ARRAYMATH_FFTFACTORY_H

#include "FFT.h"
#include "common/CPUFeatureDetector.h"
#include "common/Types.h"

namespace arraymath {

/// \brief A class for creating FFT objects.
///
/// ### Using the FFTFactory interface
///
/// \code
///    #include <ArrayMath.h>
///    #include <FFTFactory.h>
///
///    void myFunction {
///      // Initialize the FFTFactory object.
///      arraymath::FFTFactory fftFactory;
///
///      // Create an FFT object for doing 256-sized transforms.
///      const unsigned len = 256;
///      arraymath::FFT* fft = fftFactory.createFFT(len);
///      if (!fft) {
///        std::cout << "Unable to create FFT object." << std::endl;
///        return;
///      }
///
///      // Create an input array with a sine tone.
///      ArrayMath math;
///      float x[len];
///      math.ramp(x, 0.0f, 30.0f, len);
///      math.sin(x, x, len);
///
///      // Calculate the forward transform of the input signal.
///      float yReal[len], yImag[len];
///      fft->forward(yReal, yImag, x, len);
///
///      // Delete the filter object.
///      delete fft;
///    }
/// \endcode
class FFTFactory {
 public:
  /// Create a new FFT object.
  /// \param[in] size The size of the transform.
  /// \return A new FFT object, or NULL if none could be created.
  FFT* createFFT(size_t size);

 private:
  CPUFeatureDetector m_cpu;
};

}  // namespace arraymath

#endif // _ARRAYMATH_FFTFACTORY_H
