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

#include "test/TestArrayMath.h"
#include "test/TestFilterFactory.h"
#include "test/TestFFTFactory.h"
#include "test/TestRandom.h"

int main() {
  std::cout << "Testing ArrayMath..." << std::endl;
  test::testArrayMath();

  std::cout << std::endl << "Testing ArrayMath: Random..." << std::endl;
  test::testRandom();

  std::cout << std::endl << "Testing FilterFactory..." << std::endl;
  test::testFilterFactory();

  std::cout << std::endl << "Testing FFTFactory..." << std::endl;
  test::testFFTFactory();

  return 0;
}
