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
#include "common/Architecture.h"

#ifdef AM_OS_WINDOWS
# include <windows.h>
#else
# include <stdint.h>
# include <sys/time.h>
#endif

namespace arraymath {

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

const size_t kArrayLengths[] = {
  16, 128, 1024, 4096, 16384, 65536, 1048576
};

const int kNumArrayLengths = sizeof(kArrayLengths) / sizeof(size_t);

template <class OP>
void benchmark_sa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> y(arrayLength);
  math.ramp(&y[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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

void benchmark_madd_saa(ArrayMath& math, size_t arrayLength) {
  Timer timer;

  std::vector<float32> dst(arrayLength);
  std::vector<float32> y(arrayLength);
  std::vector<float32> z(arrayLength);
  math.ramp(&y[0], 1.0f, 2.0f, arrayLength);
  math.ramp(&z[0], 1.0f, 4.0f, arrayLength);

  size_t totalSamples = 0;
  double totalT = 0.0, maxSpeed = 0.0;
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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
  while (totalT < 0.2) {
    double t0 = timer.GetTime();
    size_t samples = 0;
    while (samples < 500000) {
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

#if 0
// TODO:
float32 (*p_max_f32)(const float32 *x, size_t length);
float32 (*p_min_f32)(const float32 *x, size_t length);
void (*p_pow_f32_as)(float32 *dst, const float32 *x, float32 y, size_t length);
void (*p_pow_f32_aa)(float32 *dst, const float32 *x, const float32 *y, size_t length);
void (*p_round_f32)(float32 *dst, const float32 *x, size_t length);
void (*p_sin_f32)(float32 *dst, const float32 *x, size_t length);
void (*p_sqrt_f32)(float32 *dst, const float32 *x, size_t length);
void (*p_tan_f32)(float32 *dst, const float32 *x, size_t length);
void (*p_clamp_f32)(float32 *dst, const float32 *x, float32 xMin, float32 xMax, size_t length);
void (*p_fract_f32)(float32 *dst, const float32 *x, size_t length);
void (*p_ramp_f32)(float32 *dst, float32 first, float32 last, size_t length);
void (*p_sign_f32)(float32 *dst, const float32 *x, size_t length);
float32 (*p_sum_f32)(const float32 *x, size_t length);
void (*p_sampleLinear_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
void (*p_sampleLinearRepeat_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
void (*p_sampleCubic_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
void (*p_sampleCubicRepeat_f32)(float32 *dst, const float32 *x, const float32 *t, size_t length, size_t xLength);
#endif

}

} // namespace arraymath

int main() {
  arraymath::benchmarkArrayMath();

  return 0;
}
