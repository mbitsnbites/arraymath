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

#define TEST_A(NAME, XVALUES, RESULTS, METHOD, COMPARE) \
  do { \
    beginTest(NAME); \
    AlignedArray dst(kArraySize), x(kArraySize); \
    static const size_t numValues = sizeof(RESULTS) / sizeof(RESULTS[0]); \
    for (size_t size = 1; size <= kMaxArraySize; ++size) { \
      for (size_t j = 0; j <= kMaxUnalignment; ++j) { \
        float* dst_ = dst.get() + j; \
        float* x_ = x.get() + j; \
        for (size_t k = 0; k < numValues; ++k) { \
          fillArray(x_, size, XVALUES[k]); \
          METHOD(dst_, x_, size); \
          expectAll(dst_, size, RESULTS[k], COMPARE); \
        } \
      } \
    } \
    endTest(); \
  } while(0)

#define TEST_SA(NAME, XVALUES, YVALUES, RESULTS, METHOD, COMPARE) \
  do { \
    beginTest(NAME); \
    AlignedArray dst(kArraySize), x(kArraySize), y(kArraySize); \
    static const size_t numValues = sizeof(RESULTS) / sizeof(float); \
    for (size_t size = 1; size <= kMaxArraySize; ++size) { \
      for (size_t j = 0; j <= kMaxUnalignment; ++j) { \
        float* dst_ = dst.get() + j; \
        float* y_ = y.get() + j; \
        for (size_t k = 0; k < numValues; ++k) { \
          fillArray(y_, size, YVALUES[k]); \
          METHOD(dst_, XVALUES[k], y_, size); \
          expectAll(dst_, size, RESULTS[k], COMPARE); \
        } \
      } \
    } \
    endTest(); \
  } while(0)

#define TEST_AA(NAME, XVALUES, YVALUES, RESULTS, METHOD, COMPARE) \
  do { \
    beginTest(NAME); \
    AlignedArray dst(kArraySize), x(kArraySize), y(kArraySize); \
    static const size_t numValues = sizeof(RESULTS) / sizeof(float); \
    for (size_t size = 1; size <= kMaxArraySize; ++size) { \
      for (size_t j = 0; j <= kMaxUnalignment; ++j) { \
        float* dst_ = dst.get() + j; \
        float* x_ = x.get() + j; \
        float* y_ = y.get() + j; \
        for (size_t k = 0; k < numValues; ++k) { \
          fillArray(x_, size, XVALUES[k]); \
          fillArray(y_, size, YVALUES[k]); \
          METHOD(dst_, x_, y_, size); \
          expectAll(dst_, size, RESULTS[k], COMPARE); \
        } \
      } \
    } \
    endTest(); \
  } while(0)

#define CALC_RESULTS_1(RESULTS, XVALUES, OP) \
  const size_t len_ = sizeof(XVALUES) / sizeof(XVALUES[0]); \
  float RESULTS[len_]; \
  for (size_t i_ = 0; i_ < len_; ++i_) \
    RESULTS[i_] = OP(XVALUES[i_]);

class ArrayMathTester : public Tester {
  private:
    static const size_t kMaxUnalignment = 15;
    static const size_t kMaxArraySize = 150;
    static const size_t kArraySize = kMaxArraySize + kMaxUnalignment;

    arraymath::ArrayMath m_math;

    void fillArray(float* array, size_t size, float value) {
      for (size_t i = 0; i < size; ++i)
        array[i] = value;
    }

  public:
    void runTests() {
      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, 5.0f, 5.0f };
        static const float results[] = { 5.0f, 6.0f, 1001015.0f };
        TEST_SA("add - scalar", xValues, yValues, results, m_math.add, compareExact);
        TEST_AA("add - vector", xValues, yValues, results, m_math.add, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, 5.0f, 5.0f };
        static const float results[] = { -5.0f, -4.0f, 1001005.0f };
        TEST_SA("sub - scalar", xValues, yValues, results, m_math.sub, compareExact);
        TEST_AA("sub - vector", xValues, yValues, results, m_math.sub, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 1.0f, 1001010.0f };
        static const float yValues[] = { 5.0f, -5.0f, 5.0f };
        static const float results[] = { 0.0f, -5.0f, 5005050.0f };
        TEST_SA("mul - scalar", xValues, yValues, results, m_math.mul, compareExact);
        TEST_AA("mul - vector", xValues, yValues, results, m_math.mul, compareExact);
      }

      {
        // TODO(m): mulCplx
      }

      {
        static const float xValues[] = { 0.0f, -5.0f, 1001010.0f, 4.0f };
        static const float yValues[] = { 5.0f, 2.0f, 7.0f, 3.0f };
        static const float results[] = { 0.0f, -2.5f, 143001.429f, 1.3333334f };
        TEST_SA("div - scalar", xValues, yValues, results, m_math.div, compare23bit);
        TEST_AA("div - vector", xValues, yValues, results, m_math.div, compare23bit);
      }

      {
        // TODO(m): divCplx
      }

      {
        // TODO(m): madd
      }

      {
        static const float xValues[] = { 0.0f, -5.0f, 1001010.0f };
        CALC_RESULTS_1(results, xValues, std::abs);
        TEST_A("abs", xValues, results, m_math.abs, compareExact);
      }

      {
        // TODO(m): absCplx
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        CALC_RESULTS_1(results, xValues, std::acos);
        TEST_A("acos", xValues, results, m_math.acos, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        CALC_RESULTS_1(results, xValues, std::asin);
        TEST_A("asin", xValues, results, m_math.asin, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 1001010.0f };
        CALC_RESULTS_1(results, xValues, std::atan);
        TEST_A("atan", xValues, results, m_math.atan, compareExact);
      }

      {
        // TODO(m): atan2
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        CALC_RESULTS_1(results, xValues, std::ceil);
        TEST_A("ceil", xValues, results, m_math.ceil, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        CALC_RESULTS_1(results, xValues, std::cos);
        TEST_A("cos", xValues, results, m_math.cos, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        CALC_RESULTS_1(results, xValues, std::exp);
        TEST_A("exp", xValues, results, m_math.exp, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        CALC_RESULTS_1(results, xValues, std::floor);
        TEST_A("floor", xValues, results, m_math.floor, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        CALC_RESULTS_1(results, xValues, std::log);
        TEST_A("log", xValues, results, m_math.log, compare23bit);
      }

      {
        // TODO(m): max
      }

      {
        // TODO(m): min
      }

      {
        // TODO(m): pow
      }

      {
        // TODO(m): We need to define what should happen for -2.5, -1.5, 0.5,
        // 1.5, 2.5, etc.
        static const float xValues[] = { 0.0f, 0.2f, 0.7f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const float results[] = { 0.0f, 0.0f, 1.0f, 0.0f, -1.0f, -1.0f, 10000000000.0f };
        TEST_A("round", xValues, results, m_math.round, compareExact);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        CALC_RESULTS_1(results, xValues, std::sin);
        TEST_A("sin", xValues, results, m_math.sin, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, 0.2f, 0.7f, -0.5f, 9.0f, 81.25f, 10000000000.0f };
        CALC_RESULTS_1(results, xValues, std::sqrt);
        TEST_A("sqrt", xValues, results, m_math.sqrt, compare23bit);
      }

      {
        static const float xValues[] = { 0.0f, 0.5f, -1.0f, 3.14159265f, 10000.0f, -1000.0 };
        CALC_RESULTS_1(results, xValues, std::tan);
        TEST_A("tan", xValues, results, m_math.tan, compare23bit);
      }

      {
        // TODO(m): clamp
      }

      {
        static const float xValues[] = { 0.0f, 0.2f, 3.375f, -0.25f, -7.5f, -1.0f, 10000000000.125f };
        static const float results[] = { 0.0f, 0.2f, 0.375f, 0.75f, 0.5f, 0.0f, 0.0f };
        TEST_A("fract", xValues, results, m_math.fract, compareExact);
      }

      {
        // TODO(m): fill
      }

      {
        // TODO(m): ramp
      }

      {
        static const float xValues[] = { 0.0f, 0.2f, 0.7f, -0.2f, -0.7f, -1.0f, 10000000000.0f };
        static const float results[] = { 1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, 1.0f };
        TEST_A("sign", xValues, results, m_math.sign, compareExact);
      }

      {
        // TODO(m): sum
      }

      // TODO(m): Add unit tests for random() and sample*() methods. Do it
      // in another test unit?
    }
};

void testArrayMath() {
  ArrayMathTester tester;
  tester.runTests();
}

} // namespace test
