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

#include "test/TestArrayMath.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "ArrayMath.h"
#include "test/AlignedArray.h"
#include "test/Tester.h"

namespace test {

class ArrayMathTester : public Tester {
  private:
    static const size_t kMaxArraySize = 150;  // Array length.
    static const size_t kMaxUnalignment = 15; // Unalignment of destination.
    static const size_t kMaxMisalignment = 2; // Misallignment (src vs dst).

    // The size to allocate for each test, accomodating for different alignments
    static const size_t kArraySize = kMaxArraySize + kMaxUnalignment + kMaxMisalignment;

    arraymath::ArrayMath m_math;

    void calcResults(float* results, const float* xValues, size_t numValues, float (*op)(float)) {
      for (size_t i = 0; i < numValues; ++i)
        results[i] = op(xValues[i]);
    }

    void testA(const char* name,
                const float* xValues,
                const float* results,
                size_t numValues,
                void (arraymath::ArrayMath::*method)(arraymath::float32*, const arraymath::float32*, size_t) const,
                bool (*compare)(float, float)) {
      beginTest(name);
      AlignedArray dst(kArraySize), x(kArraySize);
      for (size_t size = 1; size <= kMaxArraySize; ++size) {
        for (size_t u = 0; u <= kMaxUnalignment; ++u) {
          float* dstPtr = dst.get() + u;
          for (size_t m = 0; m <= kMaxMisalignment; ++m) {
            float* xPtr = x.get() + u + m;
            for (size_t k = 0; k < numValues; ++k) {
              fillArray(dstPtr, size, 42.0f);
              fillArray(xPtr, size, xValues[k]);
              (m_math.*method)(dstPtr, xPtr, size);
              expectAll(dstPtr, size, results[k], compare);
            }
          }
        }
      }
      endTest();
    }

    void testSA(const char* name,
                const float* xValues,
                const float* yValues,
                const float* results,
                size_t numValues,
                void (arraymath::ArrayMath::*method)(arraymath::float32*, arraymath::float32, const arraymath::float32*, size_t) const,
                bool (*compare)(float, float)) {
      beginTest(name);
      AlignedArray dst(kArraySize), y(kArraySize);
      for (size_t size = 1; size <= kMaxArraySize; ++size) {
        for (size_t u = 0; u <= kMaxUnalignment; ++u) {
          float* dstPtr = dst.get() + u;
          for (size_t m = 0; m <= kMaxMisalignment; ++m) {
            float* yPtr = y.get() + u + m;
            for (size_t k = 0; k < numValues; ++k) {
              fillArray(dstPtr, size, 42.0f);
              fillArray(yPtr, size, yValues[k]);
              (m_math.*method)(dstPtr, xValues[k], yPtr, size);
              expectAll(dstPtr, size, results[k], compare);
            }
          }
        }
      }
      endTest();
    }

    void testAA(const char* name,
                const float* xValues,
                const float* yValues,
                const float* results,
                size_t numValues,
                void (arraymath::ArrayMath::*method)(arraymath::float32*, const arraymath::float32*, const arraymath::float32*, size_t) const,
                bool (*compare)(float, float)) {
      beginTest(name);
      AlignedArray dst(kArraySize), x(kArraySize), y(kArraySize);
      for (size_t size = 1; size <= kMaxArraySize; ++size) {
        for (size_t u = 0; u <= kMaxUnalignment; ++u) {
          float* dstPtr = dst.get() + u;
          for (size_t m1 = 0; m1 <= kMaxMisalignment; ++m1) {
            float* xPtr = x.get() + u + m1;
            for (size_t m2 = 0; m2 <= kMaxMisalignment; ++m2) {
              float* yPtr = y.get() + u + m1 + m2;
              for (size_t k = 0; k < numValues; ++k) {
                fillArray(dstPtr, size, 42.0f);
                fillArray(xPtr, size, xValues[k]);
                fillArray(yPtr, size, yValues[k]);
                (m_math.*method)(dstPtr, xPtr, yPtr, size);
                expectAll(dstPtr, size, results[k], compare);
              }
            }
          }
        }
      }
      endTest();
    }

    void testSACplx(const char* name,
                    const float* xReValues,
                    const float* xImValues,
                    const float* yReValues,
                    const float* yImValues,
                    const float* resultsRe,
                    const float* resultsIm,
                    size_t numValues,
                    void (arraymath::ArrayMath::*method)(arraymath::float32*, arraymath::float32*, arraymath::float32, arraymath::float32, const arraymath::float32*, const arraymath::float32*, size_t) const,
                    bool (*compare)(float, float)) {
      beginTest(name);
      AlignedArray dstRe(kArraySize), dstIm(kArraySize), yRe(kArraySize), yIm(kArraySize);
      for (size_t size = 1; size <= kMaxArraySize; ++size) {
        for (size_t u = 0; u <= kMaxUnalignment; ++u) {
          float* dstRePtr = dstRe.get() + u;
          float* dstImPtr = dstIm.get() + u;
          for (size_t m = 0; m <= kMaxMisalignment; ++m) {
            float* yRePtr = yRe.get() + u + m;
            float* yImPtr = yIm.get() + u + m;
            for (size_t k = 0; k < numValues; ++k) {
              fillArray(dstRePtr, size, 42.0f);
              fillArray(dstImPtr, size, 42.0f);
              fillArray(yRePtr, size, yReValues[k]);
              fillArray(yImPtr, size, yImValues[k]);
              (m_math.*method)(dstRePtr, dstImPtr, xReValues[k], xImValues[k], yRePtr, yImPtr, size);
              expectAll(dstRePtr, size, resultsRe[k], compare);
              expectAll(dstImPtr, size, resultsIm[k], compare);
            }
          }
        }
      }
      endTest();
    }

    void testAACplx(const char* name,
                    const float* xReValues,
                    const float* xImValues,
                    const float* yReValues,
                    const float* yImValues,
                    const float* resultsRe,
                    const float* resultsIm,
                    size_t numValues,
                    void (arraymath::ArrayMath::*method)(arraymath::float32*, arraymath::float32*, const arraymath::float32*, const arraymath::float32*, const arraymath::float32*, const arraymath::float32*, size_t) const,
                    bool (*compare)(float, float)) {
      beginTest(name);
      AlignedArray dstRe(kArraySize), dstIm(kArraySize), xRe(kArraySize), xIm(kArraySize), yRe(kArraySize), yIm(kArraySize);
      for (size_t size = 1; size <= kMaxArraySize; ++size) {
        for (size_t u = 0; u <= kMaxUnalignment; ++u) {
          float* dstRePtr = dstRe.get() + u;
          float* dstImPtr = dstIm.get() + u;
          for (size_t m1 = 0; m1 <= kMaxMisalignment; ++m1) {
            float* xRePtr = xRe.get() + u + m1;
            float* xImPtr = xIm.get() + u + m1;
            for (size_t m2 = 0; m2 <= kMaxMisalignment; ++m2) {
              float* yRePtr = yRe.get() + u + m1 + m2;
              float* yImPtr = yIm.get() + u + m1 + m2;
              for (size_t k = 0; k < numValues; ++k) {
                fillArray(dstRePtr, size, 42.0f);
                fillArray(dstImPtr, size, 42.0f);
                fillArray(xRePtr, size, xReValues[k]);
                fillArray(xImPtr, size, xImValues[k]);
                fillArray(yRePtr, size, yReValues[k]);
                fillArray(yImPtr, size, yImValues[k]);
                (m_math.*method)(dstRePtr, dstImPtr, xRePtr, xImPtr, yRePtr, yImPtr, size);
                expectAll(dstRePtr, size, resultsRe[k], compare);
                expectAll(dstImPtr, size, resultsIm[k], compare);
              }
            }
          }
        }
      }
      endTest();
    }


  public:
    void runTests() {
      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, 5.0f, 5.0f };
        static const float results[] = { 5.0f, 6.0f, 1001015.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testSA("add - scalar", xValues, yValues, results, numValues, &arraymath::ArrayMath::add, compareExact);
        testAA("add - vector", xValues, yValues, results, numValues, &arraymath::ArrayMath::add, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, 5.0f, 5.0f };
        static const float results[] = { -5.0f, -4.0f, 1001005.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testSA("sub - scalar", xValues, yValues, results, numValues, &arraymath::ArrayMath::sub, compareExact);
        testAA("sub - vector", xValues, yValues, results, numValues, &arraymath::ArrayMath::sub, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, -5.0f, 5.0f };
        static const float results[] = { 0.0f, -5.0f, 5005050.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testSA("mul - scalar", xValues, yValues, results, numValues, &arraymath::ArrayMath::mul, compareExact);
        testAA("mul - vector", xValues, yValues, results, numValues, &arraymath::ArrayMath::mul, compareExact);
      }

      {
        static const float xReValues[] = { 0.0f, 1.0f, 0.0f, 0.0f };
        static const float xImValues[] = { 0.0f, 0.0f, 2.0f, 1001010.0f };
        static const float yReValues[] = { 5.0f, -5.0f, 3.0f, 0.0f };
        static const float yImValues[] = { 8.0f, 0.0f, 0.0f, 5.0f };
        static const float resultsRe[] = { 0.0f, -5.0f, 0.0f, -5005050.0f };
        static const float resultsIm[] = { 0.0f, 0.0f, 6.0f, 0.0f };
        static const size_t numValues = sizeof(xReValues) / sizeof(xReValues[0]);
        testSACplx("mulCplx - scalar", xReValues, xImValues, yReValues, yImValues, resultsRe, resultsIm, numValues, &arraymath::ArrayMath::mulCplx, compareExact);
        testAACplx("mulCplx - vector", xReValues, xImValues, yReValues, yImValues, resultsRe, resultsIm, numValues, &arraymath::ArrayMath::mulCplx, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, -5.0f, 1001010.0f, 4.0f };
        static const float yValues[] = { 5.0f, 2.0f, 7.0f, 3.0f };
        static const float results[] = { 0.0f, -2.5f, 143001.429f, 1.3333334f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testSA("div - scalar", xValues, yValues, results, numValues, &arraymath::ArrayMath::div, compare23bit);
        testAA("div - vector", xValues, yValues, results, numValues, &arraymath::ArrayMath::div, compare23bit);
      }

      {
        static const float xReValues[] = { 0.0f, 5.0f, 0.0f, -5005050.0f };
        static const float xImValues[] = { 0.0f, 0.0f, 4.0f, 0.0f };
        static const float yReValues[] = { 5.0f, -1.0f, 2.0f, 0.0f };
        static const float yImValues[] = { 8.0f, 0.0f, 0.0f, 5.0f };
        static const float resultsRe[] = { 0.0f, -5.0f, 0.0f, 0.0f };
        static const float resultsIm[] = { 0.0f, 0.0f, 2.0f, 1001010.0f };
        static const size_t numValues = sizeof(xReValues) / sizeof(xReValues[0]);
        testSACplx("divCplx - scalar", xReValues, xImValues, yReValues, yImValues, resultsRe, resultsIm, numValues, &arraymath::ArrayMath::divCplx, compare23bit);
        testAACplx("divCplx - vector", xReValues, xImValues, yReValues, yImValues, resultsRe, resultsIm, numValues, &arraymath::ArrayMath::divCplx, compare23bit);
      }

      {
        // TODO(m): madd
      }

      {
        static const float xValues[] = { 0.0f, -5.0f, 1001010.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::abs);
        testA("abs", xValues, results, numValues, &arraymath::ArrayMath::abs, compareExact);
      }

      {
        // TODO(m): absCplx
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::acos);
        testA("acos", xValues, results, numValues, &arraymath::ArrayMath::acos, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::asin);
        testA("asin", xValues, results, numValues, &arraymath::ArrayMath::asin, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::atan);
        testA("atan", xValues, results, numValues, &arraymath::ArrayMath::atan, compareExact);
      }

      {
        // TODO(m): atan2
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::ceil);
        testA("ceil", xValues, results, numValues, &arraymath::ArrayMath::ceil, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::cos);
        testA("cos", xValues, results, numValues, &arraymath::ArrayMath::cos, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::exp);
        testA("exp", xValues, results, numValues, &arraymath::ArrayMath::exp, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::floor);
        testA("floor", xValues, results, numValues, &arraymath::ArrayMath::floor, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::log);
        testA("log", xValues, results, numValues, &arraymath::ArrayMath::log, compare23bit);
      }

      {
        // TODO(m): max
      }

      {
        // TODO(m): min
      }

      {
        static const float xValues[] = { 999.0f, 2.2f, 7.375f, 0.25f, -1.0f };
        static const float yValues[] = { 0.0f, -1.0f, 0.5f, -1.25f, -1.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        for (size_t i = 0; i < numValues; ++i)
          results[i] = std::pow(xValues[i], yValues[i]);
        {
          beginTest("pow - scalar");
          AlignedArray dst(kArraySize), x(kArraySize);
          for (size_t size = 1; size <= kMaxArraySize; ++size) {
            for (size_t u = 0; u <= kMaxUnalignment; ++u) {
              float* dstPtr = dst.get() + u;
              for (size_t m = 0; m <= kMaxMisalignment; ++m) {
                float* xPtr = x.get() + u + m;
                for (size_t k = 0; k < sizeof(xValues) / sizeof(xValues[0]); ++k) {
                  fillArray(dstPtr, size, 42.0f);
                  fillArray(xPtr, size, xValues[k]);
                  m_math.pow(dstPtr, xPtr, yValues[k], size);
                  expectAll(dstPtr, size, results[k], compare23bit);
                }
              }
            }
          }
          endTest();
        }
        {
          beginTest("pow - array");
          AlignedArray dst(kArraySize), x(kArraySize), y(kArraySize);
          for (size_t size = 1; size <= kMaxArraySize; ++size) {
            for (size_t u = 0; u <= kMaxUnalignment; ++u) {
              float* dstPtr = dst.get() + u;
              for (size_t m1 = 0; m1 <= kMaxMisalignment; ++m1) {
                float* xPtr = x.get() + u + m1;
                for (size_t m2 = 0; m2 <= kMaxMisalignment; ++m2) {
                  float* yPtr = y.get() + u + m1 + m2;
                  for (size_t k = 0; k < sizeof(xValues) / sizeof(xValues[0]); ++k) {
                    fillArray(dstPtr, size, 42.0f);
                    fillArray(xPtr, size, xValues[k]);
                    fillArray(yPtr, size, yValues[k]);
                    m_math.pow(dstPtr, xPtr, yPtr, size);
                    expectAll(dstPtr, size, results[k], compare23bit);
                  }
                }
              }
            }
          }
          endTest();
        }
      }

      {
        // TODO(m): We need to define what should happen for -2.5, -1.5, 0.5,
        // 1.5, 2.5, etc.
        static const float xValues[] = { 0.0f, 0.2f, 0.7f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const float results[] = { 0.0f, 0.0f, 1.0f, 0.0f, -1.0f, -1.0f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testA("round", xValues, results, numValues, &arraymath::ArrayMath::round, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::sin);
        testA("sin", xValues, results, numValues, &arraymath::ArrayMath::sin, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::sqrt);
        testA("sqrt", xValues, results, numValues, &arraymath::ArrayMath::sqrt, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        float results[numValues];
        calcResults(results, xValues, numValues, std::tan);
        testA("tan", xValues, results, numValues, &arraymath::ArrayMath::tan, compare23bit);
      }

      {
        static const float xValues[] = { 1.5f, -2.2f, 7.375f, -0.25f, -1.0f, 20000000000.0f };
        static const float minValues[] = { 0.0f, -1.0f, -3.375f, -1.25f, -1.0f, 0.0f };
        static const float maxValues[] = { 4.0f, 1.0f, 3.375f, -0.5f, -1.0f, 10000000000.0f };
        static const float results[] = { 1.5f, -1.0f, 3.375f, -0.5f, -1.0f, 10000000000.0f };
        beginTest("clamp");
        AlignedArray dst(kArraySize), x(kArraySize);
        for (size_t size = 1; size <= kMaxArraySize; ++size) {
          for (size_t u = 0; u <= kMaxUnalignment; ++u) {
            float* dstPtr = dst.get() + u;
            for (size_t m = 0; m <= kMaxMisalignment; ++m) {
              float* xPtr = x.get() + u + m;
              for (size_t k = 0; k < sizeof(xValues) / sizeof(xValues[0]); ++k) {
                fillArray(dstPtr, size, 42.0f);
                fillArray(xPtr, size, xValues[k]);
                m_math.clamp(dstPtr, xPtr, minValues[k], maxValues[k], size);
                expectAll(dstPtr, size, results[k], compareExact);
              }
            }
          }
        }
        endTest();
      }

      {
        static const float xValues[] = { 0.0f, 0.2f, 3.375f, -0.25f, -7.5f, -1.0f, 10000000000.125f };
        static const float results[] = { 0.0f, 0.2f, 0.375f, 0.75f, 0.5f, 0.0f, 0.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testA("fract", xValues, results, numValues, &arraymath::ArrayMath::fract, compareExact);
      }

      {
        static const float values[] = { 0.0f, 0.2f, 3.375f, -0.25f, -7.5f, -1.0f, 10000000000.125f };
        beginTest("fill");
        AlignedArray dst(kArraySize);
        for (size_t size = 1; size <= kMaxArraySize; ++size) {
          for (size_t u = 0; u <= kMaxUnalignment; ++u) {
            float* dstPtr = dst.get() + u;
            for (size_t k = 0; k < sizeof(values) / sizeof(values[0]); ++k) {
              fillArray(dstPtr, size, 42.0f);
              m_math.fill(dstPtr, values[k], size);
              expectAll(dstPtr, size, values[k], compareExact);
            }
          }
        }
        endTest();
      }

      {
        // TODO(m): ramp
      }

      {
        static const float xValues[] = { 0.0f, 0.2f, 0.7f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const float results[] = { 1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, 1.0f };
        static const size_t numValues = sizeof(xValues) / sizeof(xValues[0]);
        testA("sign", xValues, results, numValues, &arraymath::ArrayMath::sign, compareExact);
      }

      {
        // TODO(m): sum
      }

      // TODO(m): Add unit tests for sample*() methods. Do it
      // in another test unit?
    }
};

void testArrayMath() {
  ArrayMathTester tester;
  tester.runTests();
}

} // namespace test
