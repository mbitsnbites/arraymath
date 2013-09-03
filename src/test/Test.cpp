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

void testArrayMath() {
  std::cout << std::endl << "Testing ArrayMath..." << std::endl;

  arraymath::ArrayMath math;

  std::vector<float> x(100000), y(10000), z(123123);

  math.ramp(&x[0], 0.0f, 1.0f, x.size());
  math.mul(&x[0], 3.141592654f, &x[0], x.size());
  math.sin(&x[0], &x[0], x.size());

  math.ramp(&y[0], 2.0f, 8.0f, y.size());
  math.madd(&y[0], 0.5f, &x[0], &y[0], std::min(x.size(), y.size()));
  math.add(&y[0], &x[0], &y[0], std::min(x.size(), y.size()));
  math.div(&y[0], 3.5f, &y[0], y.size());

  math.random(&z[0], 0.0f, 10.0f, z.size());

  std::cout << "sum(x) = " << math.sum(&x[0], x.size()) << std::endl;
  std::cout << "min(x) = " << math.min(&x[0], x.size()) << std::endl;
  std::cout << "max(x) = " << math.max(&x[0], x.size()) << std::endl;

  std::cout << "sum(y) = " << math.sum(&y[0], y.size()) << std::endl;
  std::cout << "min(y) = " << math.min(&y[0], y.size()) << std::endl;
  std::cout << "max(y) = " << math.max(&y[0], y.size()) << std::endl;

  std::cout << "mean(z) = " << (math.sum(&z[0], z.size()) / static_cast<float>(z.size())) << std::endl;
  std::cout << "min(z) = " << math.min(&z[0], z.size()) << std::endl;
  std::cout << "max(z) = " << math.max(&z[0], z.size()) << std::endl;

  std::cout << std::endl << "sqrt([2.0, 4.0])" << std::endl;
  math.ramp(&x[0], 2.0f, 4.0f, x.size());
  math.sqrt(&x[0], &x[0], x.size());
  std::cout << "min(x) = " << math.min(&x[0], x.size()) << std::endl;
  std::cout << "max(x) = " << math.max(&x[0], x.size()) << std::endl;
}

void testFilterFactory() {
  std::cout << std::endl << "Testing FilterFactory..." << std::endl;

  arraymath::ArrayMath math;
  arraymath::FilterFactory filterFactory;

  arraymath::Filter* f = filterFactory.createFilter(8, 0);
  if (!f) {
    std::cout << "Unable to create filter." << std::endl;
    return;
  }

  arraymath::float32 b[8] = { 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f };
  f->setB(b, 8);

  std::vector<float> x(100000), y(100000);
  math.random(&x[0], -1.0f, 1.0f, x.size());
  f->filter(&y[0], &x[0], x.size());

  std::cout << "sum(y) = " << math.sum(&y[0], y.size()) << std::endl;
  std::cout << "min(y) = " << math.min(&y[0], y.size()) << std::endl;
  std::cout << "max(y) = " << math.max(&y[0], y.size()) << std::endl;

  delete f;
}

int main() {
  testArrayMath();
  testFilterFactory();

  return 0;
}
