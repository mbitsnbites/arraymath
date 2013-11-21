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

#include "test/TestFilterFactory.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "FilterFactory.h"
#include "test/AlignedArray.h"
#include "test/ScopedPtr.h"
#include "test/Tester.h"

namespace test {

class FilterFactoryTester : public Tester {
  private:
    static const size_t kMaxArraySize = 150; // Array length.
    static const int kMaxFIRSize = 50;       // FIR-filter length.
    static const int kMaxIIRSize = 8;        // IIR-filter length.

    arraymath::FilterFactory m_factory;

  public:
    void runTests() {
      {
        beginTest("FIR");
        AlignedArray dst(kMaxArraySize), x(kMaxArraySize);
        for (int size = 1; size <= kMaxFIRSize; ++size) {
          for (size_t length = 1; length <= kMaxArraySize; ++length) {
            ScopedPtr<arraymath::Filter> f(m_factory.createFilter(size, 0));
            size_t passes = 2 + size / length;
            for (size_t pass = 1; pass <= passes; ++pass) {
              float value = static_cast<float>(pass + length);
              fillArray(dst.get(), length, 42.0f);
              fillArray(x.get(), length, value);
              f->filter(dst.get(), x.get(), length);
              expectAll(dst.get(), length, value, compareExact);
            }
          }
        }
        endTest();
      }

      {
        beginTest("IIR");
        AlignedArray dst(kMaxArraySize), x(kMaxArraySize);
        for (int size = 1; size <= kMaxIIRSize; ++size) {
          for (size_t length = 1; length <= kMaxArraySize; ++length) {
            ScopedPtr<arraymath::Filter> f(m_factory.createFilter(size, size - 1));
            size_t passes = 2 + size / length;
            for (size_t pass = 1; pass <= passes; ++pass) {
              float value = static_cast<float>(pass + length);
              fillArray(dst.get(), length, 42.0f);
              fillArray(x.get(), length, value);
              f->filter(dst.get(), x.get(), length);
              expectAll(dst.get(), length, value, compareExact);
            }
          }
        }
        endTest();
      }
    }
};

void testFilterFactory() {
  FilterFactoryTester tester;
  tester.runTests();
}

} // namespace test
