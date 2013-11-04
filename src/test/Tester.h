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
      return a == b;
    }

    static bool compare23bit(float a, float b) {
      float err = std::abs(a - b);
      if (err == 0.0f)
        return true;
      float denom = std::abs(b);
      if (denom == 0.0f)
        return false;
      // 0x3F800001 - 0x3F800000 = 0.0000001192092896
      return (err / denom) < 1.2e-7;
    }

    static bool compare21bit(float a, float b) {
      float err = std::abs(a - b);
      if (err == 0.0f)
        return true;
      float denom = std::abs(b);
      if (denom == 0.0f)
        return false;
      // 0x3F800004 - 0x3F800000 = 0.0000004768371582
      return (err / denom) < 4.8e-7;
    }

  private:
    int m_failCount;
    int m_passCount;
};

} // namespace test

#endif // _ARRAYMATH_TEST_TESTER_H
