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

#ifndef _ARRAYMATH_TEST_TESTER_H
#define _ARRAYMATH_TEST_TESTER_H

#include <cmath>
#include <cstddef>

namespace test {

class Tester {
  protected:
    void beginTest(const char* name);
    void endTest();

    bool expectAll(const float* array, size_t size, float value, bool (*compare)(float, float));

    static bool compareExact(float a, float b) {
      return a == b || (std::isnan(a) && std::isnan(b));
    }

    // Mantissa precision factors:
    // 23-bit: 0x3f800001 - 0x3f800000 = 0.0000001192092896
    // 22-bit: 0x3f800002 - 0x3f800000 = 0.000000238418579
    // 21-bit: 0x3f800004 - 0x3f800000 = 0.0000004768371582
    // 20-bit: 0x3f800008 - 0x3f800000 = 0.0000009536743164
    // 19-bit: 0x3f800010 - 0x3f800000 = 0.0000019073486328
    // 18-bit: 0x3f800020 - 0x3f800000 = 0.0000038146972656
    // 17-bit: 0x3f800040 - 0x3f800000 = 0.0000076293945312
    // 16-bit: 0x3f800080 - 0x3f800000 = 0.0000152587890625
    // 15-bit: 0x3f800100 - 0x3f800000 = 0.000030517578125
    // 14-bit: 0x3f800200 - 0x3f800000 = 0.00006103515625

    static bool compare23bit(float a, float b) {
      return compareWithThreshold(a, b, 1.2e-7);
    }

    static bool compare22bit(float a, float b) {
      return compareWithThreshold(a, b, 2.4e-7);
    }

    static bool compare21bit(float a, float b) {
      return compareWithThreshold(a, b, 4.8e-7);
    }

    static bool compareLT(float a, float b) {
      return a < b && !std::isnan(b);
    }

    static bool compareGE(float a, float b) {
      return a >= b && !std::isnan(b);
    }

  private:
    static bool compareWithThreshold(float a, float b, float threshold) {
      if (compareExact(a, b))
        return true;
      float denom = std::abs(b);
      if (denom == 0.0f)
        return false;
      return (std::abs(a - b) / denom) < threshold;
    }

    int m_failCount;
    int m_passCount;
};

} // namespace test

#endif // _ARRAYMATH_TEST_TESTER_H
