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
        static const float xValues[] = { 0.0f, -5.0f, 1001010.0f, 4.0f };
        static const float yValues[] = { 5.0f, 2.0f, 7.0f, 3.0f };
        static const float results[] = { 0.0f, -2.5f, 143001.429f, 1.3333334f };
        TEST_SA("div - scalar", xValues, yValues, results, m_math.div, compare23bit);
        TEST_AA("div - vector", xValues, yValues, results, m_math.div, compare23bit);
      }

      // TODO(m): Implement all unit tests.
    }
};

void testArrayMath() {
  ArrayMathTester tester;
  tester.runTests();
}

} // namespace test
