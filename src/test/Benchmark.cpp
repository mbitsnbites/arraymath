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

#include <iostream>
#include <vector>

#include "ArrayMath.h"
#include "FilterFactory.h"
#include "FFTFactory.h"
#include "common/Architecture.h"

#ifdef AM_OS_WINDOWS
# include <windows.h>
#else
# include <stdint.h>
# include <sys/time.h>
#endif

namespace arraymath {

//------------------------------------------------------------------------------
// Benchmark configuration.
//------------------------------------------------------------------------------

const size_t kArrayLengths[] = {
  16, 128, 1024, 4096, 16384, 65536, 1048576
};

const int kNumArrayLengths = sizeof(kArrayLengths) / sizeof(size_t);

const double kMinTimePerTest = 0.2;

const size_t kMinSamplesPerTest = 500000;


//------------------------------------------------------------------------------
// The timer class is used for measuring time.
//------------------------------------------------------------------------------

class Timer {
 public:
  Timer() {
#ifdef AM_OS_WINDOWS
    if(QueryPerformanceFrequency((LARGE_INTEGER *)&m_timeFreq)) {
      QueryPerformanceCounter((LARGE_INTEGER *)&m_timeStart);
    }
    else {
      m_timeFreq = 0;
    }
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    m_timeStart = 1000000 * static_cast<int64_t>(tv.tv_sec) +
                  static_cast<int64_t>(tv.tv_usec);
#endif
  }

  double GetTime() {
#ifdef AM_OS_WINDOWS
    __int64 t;
    QueryPerformanceCounter((LARGE_INTEGER *)&t);
    return double(t - m_timeStart) / double(m_timeFreq);
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    int64_t t = 1000000 * static_cast<int64_t>(tv.tv_sec) +
                static_cast<int64_t>(tv.tv_usec);
    return (1e-6) * double(t - m_timeStart);
#endif
  }

 private:
#ifdef AM_OS_WINDOWS
  __int64 m_timeFreq;
  __int64 m_timeStart;
#else
  int64_t m_timeStart;
#endif
};


//------------------------------------------------------------------------------
// Benchmark routines.
//------------------------------------------------------------------------------

template <class OP>
void benchmark_sa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> y(arrayLength);
  math.ramp(&y[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dst[0], 0.5f, &y[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_aa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  std::vector<float32> y(arrayLength);
  math.ramp(&x[0], 1.0f, 2.0f, arrayLength);
  math.ramp(&y[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dst[0], &x[0], &y[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_cplx_sa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dstReal(arrayLength), dstImag(arrayLength);
  std::vector<float32> yReal(arrayLength), yImag(arrayLength);
  math.ramp(&yReal[0], 1.0f, 4.0f, arrayLength);
  math.ramp(&yImag[0], 2.0f, 8.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dstReal[0], &dstImag[0], 0.5f, -1.0f, &yReal[0], &yImag[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_cplx_aa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dstReal(arrayLength), dstImag(arrayLength);
  std::vector<float32> xReal(arrayLength), xImag(arrayLength);
  std::vector<float32> yReal(arrayLength), yImag(arrayLength);
  math.ramp(&xReal[0], 0.5f, 2.0f, arrayLength);
  math.ramp(&xImag[0], 1.0f, 3.0f, arrayLength);
  math.ramp(&yReal[0], 1.0f, 4.0f, arrayLength);
  math.ramp(&yImag[0], 2.0f, 8.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dstReal[0], &dstImag[0], &xReal[0], &xImag[0], &yReal[0], &yImag[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_a(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  math.ramp(&x[0], 0.5f, 2.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dst[0], &x[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_const(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> x(arrayLength);
  math.ramp(&x[0], 0.5f, 2.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &x[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_madd_saa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> y(arrayLength);
  std::vector<float32> z(arrayLength);
  math.ramp(&y[0], 1.0f, 2.0f, arrayLength);
  math.ramp(&z[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.madd(&dst[0], 0.5f, &y[0], &z[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_madd_aaa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  std::vector<float32> y(arrayLength);
  std::vector<float32> z(arrayLength);
  math.ramp(&x[0], 0.5f, 1.5f, arrayLength);
  math.ramp(&y[0], 1.0f, 2.0f, arrayLength);
  math.ramp(&z[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.madd(&dst[0], &x[0], &y[0], &z[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_pow_as(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  math.ramp(&x[0], 0.5f, 2.0f, arrayLength);
  math.sin(&x[0], &x[0], arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.pow(&dst[0], &x[0], 5.0f, arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_random(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.random(&dst[0], -0.9f, 0.9f, arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_clamp(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  math.ramp(&x[0], 0.0f, 20.0f, arrayLength);
  math.sin(&x[0], &x[0], arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.clamp(&dst[0], &x[0], -0.9f, 0.9f, arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

void benchmark_ramp(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      math.ramp(&dst[0], -0.9f, 0.9f, arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

template <class OP>
void benchmark_sampler(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  size_t xLength = arrayLength + 100;
  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(xLength);
  std::vector<float32> t(arrayLength);
  math.ramp(&x[0], 1.0f, 4.0f, xLength);
  math.ramp(&t[0], -2.5f, static_cast<float32>(xLength + 2), arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      OP::op(math, &dst[0], &x[0], &t[0], arrayLength, xLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;
}

struct add_sa_OP {
  static void op(ArrayMath& math, float32* dst, float32 x, const float32* y, unsigned length) {
    math.add(dst, x, y, length);
  }
};

struct add_aa_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.add(dst, x, y, length);
  }
};

struct sub_sa_OP {
  static void op(ArrayMath& math, float32* dst, float32 x, const float32* y, unsigned length) {
    math.sub(dst, x, y, length);
  }
};

struct sub_aa_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.sub(dst, x, y, length);
  }
};

struct mul_sa_OP {
  static void op(ArrayMath& math, float32* dst, float32 x, const float32* y, unsigned length) {
    math.mul(dst, x, y, length);
  }
};

struct mul_aa_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.mul(dst, x, y, length);
  }
};

struct mul_cplx_sa_OP {
  static void op(ArrayMath& math, float32* dstReal, float32* dstImag, float32 xReal, float32 xImag, const float32* yReal, const float32* yImag, unsigned length) {
    math.mulCplx(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }
};

struct mul_cplx_aa_OP {
  static void op(ArrayMath& math, float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag, const float32* yReal, const float32* yImag, unsigned length) {
    math.mulCplx(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }
};

struct div_sa_OP {
  static void op(ArrayMath& math, float32* dst, float32 x, const float32* y, unsigned length) {
    math.div(dst, x, y, length);
  }
};

struct div_aa_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.div(dst, x, y, length);
  }
};

struct div_cplx_sa_OP {
  static void op(ArrayMath& math, float32* dstReal, float32* dstImag, float32 xReal, float32 xImag, const float32* yReal, const float32* yImag, unsigned length) {
    math.divCplx(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }
};

struct div_cplx_aa_OP {
  static void op(ArrayMath& math, float32* dstReal, float32* dstImag, const float32* xReal, const float32* xImag, const float32* yReal, const float32* yImag, unsigned length) {
    math.divCplx(dstReal, dstImag, xReal, xImag, yReal, yImag, length);
  }
};

struct abs_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.abs(dst, x, length);
  }
};

struct absCplx_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.absCplx(dst, x, y, length);
  }
};

struct acos_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.acos(dst, x, length);
  }
};

struct asin_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.asin(dst, x, length);
  }
};

struct atan_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.atan(dst, x, length);
  }
};

struct atan2_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.atan2(dst, x, y, length);
  }
};

struct ceil_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.ceil(dst, x, length);
  }
};

struct cos_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.cos(dst, x, length);
  }
};

struct exp_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.exp(dst, x, length);
  }
};

struct floor_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.floor(dst, x, length);
  }
};

struct log_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.log(dst, x, length);
  }
};

struct pow_aa_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* y, unsigned length) {
    math.pow(dst, x, y, length);
  }
};

struct round_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.round(dst, x, length);
  }
};

struct sin_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.sin(dst, x, length);
  }
};

struct sqrt_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.sqrt(dst, x, length);
  }
};

struct tan_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.tan(dst, x, length);
  }
};

struct fract_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.fract(dst, x, length);
  }
};

struct sign_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, unsigned length) {
    math.sign(dst, x, length);
  }
};

struct max_OP {
  static float32 op(ArrayMath& math, const float32* x, unsigned length) {
    return math.max(x, length);
  }
};

struct min_OP {
  static float32 op(ArrayMath& math, const float32* x, unsigned length) {
    return math.min(x, length);
  }
};

struct sum_OP {
  static float32 op(ArrayMath& math, const float32* x, unsigned length) {
    return math.sum(x, length);
  }
};

struct sampleLinear_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* t, unsigned length, unsigned xLength) {
    math.sampleLinear(dst, x, t, length, xLength);
  }
};

struct sampleLinearRepeat_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* t, unsigned length, unsigned xLength) {
    math.sampleLinearRepeat(dst, x, t, length, xLength);
  }
};

struct sampleCubic_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* t, unsigned length, unsigned xLength) {
    math.sampleCubic(dst, x, t, length, xLength);
  }
};

struct sampleCubicRepeat_OP {
  static void op(ArrayMath& math, float32* dst, const float32* x, const float32* t, unsigned length, unsigned xLength) {
    math.sampleCubicRepeat(dst, x, t, length, xLength);
  }
};

void benchmark_filter(ArrayMath& math, FilterFactory& factory, int bSize, int aSize, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> x(arrayLength);
  math.random(&x[0], -1.0f, 1.0f, arrayLength);

  Filter* filter = factory.createFilter(bSize, aSize);
  if (!filter) {
    std::cout << "ERROR: Unable to create Filter object." << std::endl;
    return;
  }

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      filter->filter(&dst[0], &x[0], arrayLength);
      samples += arrayLength;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << arrayLength << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;

  delete filter;
}

void benchmark_fft(ArrayMath& math, FFTFactory& factory, size_t size) {
  Timer timer;

  std::vector<float32> dstReal(size);
  std::vector<float32> dstImag(size);
  std::vector<float32> x(size);
  math.random(&x[0], -1.0f, 1.0f, size);

  FFT* fft = factory.createFFT(size);
  if (!fft) {
    std::cout << "ERROR: Unable to create FFT object." << std::endl;
    return;
  }

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      fft->forward(&dstReal[0], &dstImag[0], &x[0]);
      samples += size;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << size << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;

  delete fft;
}

void benchmark_fftComplex(ArrayMath& math, FFTFactory& factory, size_t size) {
  Timer timer;

  std::vector<float32> dstReal(size);
  std::vector<float32> dstImag(size);
  std::vector<float32> xReal(size);
  std::vector<float32> xImag(size);
  math.random(&xReal[0], -1.0f, 1.0f, size);
  math.random(&xImag[0], -1.0f, 1.0f, size);

  FFT* fft = factory.createFFT(size);
  if (!fft) {
    std::cout << "ERROR: Unable to create FFT object." << std::endl;
    return;
  }

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < kMinTimePerTest) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < kMinSamplesPerTest) {
      fft->forwardCplx(&dstReal[0], &dstImag[0], &xReal[0], &xImag[0]);
      samples += size;
    }
    double dt = timer.GetTime() - t0;
    double speed = double(samples) / dt;
    if (speed > maxSpeed) {
      maxSpeed = speed;
    }
    totalT += dt;
    totalSamples += samples;
  }
  double avgSpeed = double(totalSamples) / totalT;

  std::cout << size << ": " << (1e-6 * maxSpeed) << " Msamples/s ("
            << (1e-6 * avgSpeed) << " average)" << std::endl;

  delete fft;
}


//------------------------------------------------------------------------------
// ArrayMath benhcmark.
//------------------------------------------------------------------------------

void benchmarkArrayMath() {
  std::cout << std::endl << "Benchmarking ArrayMath..." << std::endl;

  ArrayMath math;

  std::cout << std::endl << "ArrayMath.add(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sa<add_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.add()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<add_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sub(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sa<sub_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sub()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<sub_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.mul(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sa<mul_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.mul()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<mul_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.mulCplx(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_cplx_sa<mul_cplx_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.mulCplx()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_cplx_aa<mul_cplx_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.div(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sa<div_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.div()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<div_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.divCplx(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_cplx_sa<div_cplx_sa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.divCplx()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_cplx_aa<div_cplx_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.madd(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_madd_saa(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.madd()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_madd_aaa(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.abs()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<abs_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.absCplx()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<absCplx_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.acos()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<acos_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.asin()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<asin_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.atan()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<atan_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.atan2()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<atan2_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.ceil()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<ceil_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.cos()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<cos_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.exp()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<exp_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.floor()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<floor_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.log()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<log_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.max()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_const<max_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.min()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_const<min_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.pow(), scalar" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_pow_as(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.pow()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_aa<pow_aa_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.random()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_random(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.round()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<round_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sin()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<sin_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sqrt()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<sqrt_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.tan()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<tan_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.clamp()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_clamp(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.fract()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<fract_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.ramp()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_ramp(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sign()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_a<sign_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sum()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_const<sum_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sampleLinear()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sampler<sampleLinear_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sampleLinearRepeat()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sampler<sampleLinearRepeat_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sampleCubic()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sampler<sampleCubic_OP>(math, kArrayLengths[i]);
  }

  std::cout << std::endl << "ArrayMath.sampleCubicRepeat()" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_sampler<sampleCubicRepeat_OP>(math, kArrayLengths[i]);
  }
}

//------------------------------------------------------------------------------
// Filter benhcmark.
//------------------------------------------------------------------------------

void benchmarkFilter() {
  std::cout << std::endl << "Benchmarking Filter..." << std::endl;

  ArrayMath math;
  FilterFactory factory;

  std::cout << std::endl << "Filter(1,1)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 1, 1, kArrayLengths[i]);
  }

  std::cout << std::endl << "Filter(2,1)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 2, 1, kArrayLengths[i]);
  }

  std::cout << std::endl << "Filter(2,2)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 2, 2, kArrayLengths[i]);
  }

  std::cout << std::endl << "Filter(3,2)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 3, 2, kArrayLengths[i]);
  }

  std::cout << std::endl << "Filter(16,0)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 16, 0, kArrayLengths[i]);
  }

  std::cout << std::endl << "Filter(128,0)" << std::endl;
  for (int i = 0; i < kNumArrayLengths; ++i) {
    benchmark_filter(math, factory, 128, 0, kArrayLengths[i]);
  }
}

//------------------------------------------------------------------------------
// FFT benhcmark.
//------------------------------------------------------------------------------

void benchmarkFFT() {
  std::cout << std::endl << "Benchmarking FFT..." << std::endl;

  ArrayMath math;
  FFTFactory factory;

  static const size_t kFFTSizes[] = { 16, 47, 48, 49, 50, 128, 256,
                                      1000, 1024, 2048, 10000, 16384 };
  static const int kNumFFTSizes = sizeof(kFFTSizes) / sizeof(size_t);

  std::cout << std::endl << "FFT - real" << std::endl;
  for (int i = 0; i < kNumFFTSizes; ++i) {
    benchmark_fft(math, factory, kFFTSizes[i]);
  }

  std::cout << std::endl << "FFT - complex" << std::endl;
  for (int i = 0; i < kNumFFTSizes; ++i) {
    benchmark_fftComplex(math, factory, kFFTSizes[i]);
  }
}

} // namespace arraymath


int main() {
  arraymath::benchmarkArrayMath();
  arraymath::benchmarkFilter();
  arraymath::benchmarkFFT();

  return 0;
}
