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

#include <cmath>
#include <iostream>
#include <vector>

#include "ArrayMath.h"
#include "FilterFactory.h"
#include "FFTFactory.h"

namespace {

template <typename T>
void printArray(const T* v, size_t length) {
  std::cout << "[";
  if (length > 0) {
    for (size_t k = 0; k < length - 1; ++k)
      std::cout << v[k] << ", ";
    std::cout << v[length - 1];
  }
  std::cout << "]";
}

template <typename T>
void printArray(std::vector<T> v) {
  printArray(&v[0], v.size());
}

}  // anonymous namespace

void testArrayMath() {
  std::cout << std::endl << "Testing ArrayMath..." << std::endl;

  arraymath::ArrayMath math;

  std::vector<float> x(100000), y(10000), z(123123);

  math.ramp(&x[0], 0.0f, 1.0f, x.size());
  math.mul(&x[0], 3.141592654f, &x[0], x.size());
  math.sin(&x[0], &x[0], x.size());

  math.ramp(&y[0], 2.0f, 8.0f, y.size());
  math.madd(&y[0], 0.5f, &x[0], &y[0], std::min(x.size(), y.size()));
  math.add(&y[0], &x[0], &y[0], std::min(x.size(), y.size()));
  math.div(&y[0], 3.5f, &y[0], y.size());

  math.random(&z[0], 0.0f, 10.0f, z.size());

  std::cout << "sum(x) = " << math.sum(&x[0], x.size()) << std::endl;
  std::cout << "min(x) = " << math.min(&x[0], x.size()) << std::endl;
  std::cout << "max(x) = " << math.max(&x[0], x.size()) << std::endl;

  std::cout << "sum(y) = " << math.sum(&y[0], y.size()) << std::endl;
  std::cout << "min(y) = " << math.min(&y[0], y.size()) << std::endl;
  std::cout << "max(y) = " << math.max(&y[0], y.size()) << std::endl;

  std::cout << "mean(z) = " << (math.sum(&z[0], z.size()) / static_cast<float>(z.size())) << std::endl;
  std::cout << "min(z) = " << math.min(&z[0], z.size()) << std::endl;
  std::cout << "max(z) = " << math.max(&z[0], z.size()) << std::endl;

  std::cout << std::endl << "sqrt([2.0, 4.0])" << std::endl;
  math.ramp(&x[0], 2.0f, 4.0f, x.size());
  math.sqrt(&x[0], &x[0], x.size());
  std::cout << "min(x) = " << math.min(&x[0], x.size()) << std::endl;
  std::cout << "max(x) = " << math.max(&x[0], x.size()) << std::endl;

  {
    std::vector<float> x(33), y(33);
    math.ramp(&x[0], -5.0f, 5.0f, x.size());
    math.sign(&y[0], &x[0], x.size());
    std::cout << "Test sign()" << std::endl;
    std::cout << "x = "; printArray(x); std::cout << std::endl;
    std::cout << "y = "; printArray(y); std::cout << std::endl;
  }

  {
    std::vector<float> t(33), y(33), x(4);
    math.ramp(&t[0], -0.5f, 4.0f, t.size());
    math.ramp(&x[0], -1.0f, 1.0f, x.size());
    math.sampleLinear(&y[0], &x[0], &t[0], y.size(), x.size());
    std::cout << "Test sampleLinear()" << std::endl;
    std::cout << "x = "; printArray(x); std::cout << std::endl;
    std::cout << "t = "; printArray(t); std::cout << std::endl;
    std::cout << "y = "; printArray(y); std::cout << std::endl;
  }
}

void testFilterFactory() {
  std::cout << std::endl << "Testing FilterFactory..." << std::endl;

  arraymath::ArrayMath math;
  arraymath::FilterFactory filterFactory;

  arraymath::Filter* f = filterFactory.createFilter(8, 0);
  if (!f) {
    std::cout << "Unable to create filter." << std::endl;
    return;
  }

  arraymath::float32 b[8] = { 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f };
  f->setB(b, 8);

  std::vector<float> x(100000), y(100000);
  math.random(&x[0], -1.0f, 1.0f, x.size());
  f->filter(&y[0], &x[0], x.size());

  std::cout << "sum(y) = " << math.sum(&y[0], y.size()) << std::endl;
  std::cout << "min(y) = " << math.min(&y[0], y.size()) << std::endl;
  std::cout << "max(y) = " << math.max(&y[0], y.size()) << std::endl;

  delete f;
}

void testFFTFactory() {
  std::cout << std::endl << "Testing FFTFactory..." << std::endl;

  // Initialize the FFTFactory object.
  arraymath::FFTFactory fftFactory;

  // Create an FFT object for doing 256-sized transforms.
  const unsigned len = 256;
  arraymath::FFT* fft = fftFactory.createFFT(len);
  if (!fft) {
    std::cout << "Unable to create FFT object." << std::endl;
    return;
  }

  // Create an input array with a sine tone.
  arraymath::ArrayMath math;
  std::vector<float> x(len);
  math.ramp(&x[0], 0.0f, 31.41592654f, len);
  math.sin(&x[0], &x[0], len);

  // Original signal.
  {
    std::vector<float> tmp(len);
    math.mul(&tmp[0], &x[0], &x[0], len);
    float rms = std::sqrt(math.sum(&tmp[0], len) / static_cast<float>(len));
    std::cout << "rms(x) = " << rms << std::endl;
  }

  // Calculate the forward transform of the input signal.
  std::vector<float> yReal(len), yImag(len);
  fft->forward(&yReal[0], &yImag[0], &x[0]);

  // Fourier transform.
  {
    std::vector<float> tmp(len);
    math.mul(&tmp[0], &yReal[0], &yReal[0], len);
    math.madd(&tmp[0], &yImag[0], &yImag[0], &tmp[0], len);
    float rms = std::sqrt(math.sum(&tmp[0], len) / static_cast<float>(len));
    std::cout << "rms(y) = " << rms << std::endl;
  }

  // Calculate the inverse transform of the frequency signal.
  std::vector<float> x2(len);
  fft->inverse(&x2[0], &yReal[0], &yImag[0]);

  // Re-transformed signal.
  {
    std::vector<float> tmp(len);
    math.mul(&tmp[0], &x2[0], &x2[0], len);
    float rms = std::sqrt(math.sum(&tmp[0], len) / static_cast<float>(len));
    std::cout << "rms(x2) = " << rms << std::endl;
  }

  // Delete the FFT object.
  delete fft;
}

int main() {
  testArrayMath();
  testFilterFactory();
  testFFTFactory();

  return 0;
}
