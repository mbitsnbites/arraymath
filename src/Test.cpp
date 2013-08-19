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

int main() {
  arraymath::ArrayMath am;

  std::vector<float> x(100000), y(10000), z(123123);

  am.ramp(&x[0], 0.0f, 1.0f, x.size());
  am.mul(&x[0], 3.141592654f, &x[0], x.size());
  am.sin(&x[0], &x[0], x.size());

  am.ramp(&y[0], 0.0f, 1.0f, y.size());
  am.madd(&y[0], 0.5f, &x[0], &y[0], std::min(x.size(), y.size()));
  am.add(&y[0], &x[0], &y[0], std::min(x.size(), y.size()));

  am.random(&z[0], 0.0f, 10.0f, z.size());

  std::cout << "sum(x) = " << am.sum(&x[0], x.size()) << std::endl;
  std::cout << "min(x) = " << am.min(&x[0], x.size()) << std::endl;
  std::cout << "max(x) = " << am.max(&x[0], x.size()) << std::endl;

  std::cout << "sum(y) = " << am.sum(&y[0], y.size()) << std::endl;
  std::cout << "min(y) = " << am.min(&y[0], y.size()) << std::endl;
  std::cout << "max(y) = " << am.max(&y[0], y.size()) << std::endl;

  std::cout << "mean(z) = " << (am.sum(&z[0], z.size()) / static_cast<float>(z.size())) << std::endl;
  std::cout << "min(z) = " << am.min(&z[0], z.size()) << std::endl;
  std::cout << "max(z) = " << am.max(&z[0], z.size()) << std::endl;

  return 0;
}
