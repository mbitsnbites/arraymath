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

#include "test/TestRandom.h"

#include <iostream>

#include "ArrayMath.h"
#include "test/AlignedArray.h"
#include "test/Tester.h"

namespace test {

class RandomTester : public Tester {
  private:
    arraymath::ArrayMath m_math;

    static const int kNumIterations = 100;

  public:
    void runTests() {
      static const size_t kArrayLengths[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000, 100000
      };
      static const int kNumArrayLengths = sizeof(kArrayLengths) / sizeof(size_t);
      static const size_t kMaxArrayLength = kArrayLengths[kNumArrayLengths - 1];

      {
        beginTest("range");
        AlignedArray dst(kMaxArrayLength);
        for (int j = 0; j < kNumIterations; ++j) {
          for (int i = 0; i < kNumArrayLengths; ++i) {
            size_t length = kArrayLengths[i];
            m_math.random(dst.get(), 0.0f,1.0f, length);
            expectAll(dst.get(), length, 0.0f, compareGE);
            expectAll(dst.get(), length, 1.0f, compareLT);
          }
        }
        endTest();
      }
    }
};

void testRandom() {
  RandomTester tester;
  tester.runTests();
}

} // namespace test
