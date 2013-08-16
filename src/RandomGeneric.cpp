// -*- Mode: c++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//------------------------------------------------------------------------------
// ArrayMath - an array math library
//------------------------------------------------------------------------------
// Copyright (c) 2013 Marcus Geelnard
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

//------------------------------------------------------------------------------
// This is an implementation of the WELL512 pseudo random number generator. Its
// main characteristics are that it's very fast and has a very long period
// (2^512-1 samples).
//------------------------------------------------------------------------------

#include "RandomGeneric.h"

namespace arraymath {

RandomGeneric::RandomGeneric() {
  // Seed the state (16 32-bit integers).
  const unsigned seed = 5489;
  unsigned* state = m_state;
  unsigned x = *state++ = seed;
  for (unsigned i = 1; i < 16; ++i) {
    x = *state++ = 1812433253 * (x ^ (x >> 30)) + i;
  }
  m_index = 0;
}

void RandomGeneric::random(float32 *dst, float32 low, float32 high, size_t length) {
  #define MUTATE_LEFT(value, shift) value ^ (value << shift)
  #define MUTATE_RIGHT(value, shift)  value ^ (value >> shift)
  #define MUTATE_LEFT_MIX(value, shift, mix)  value ^ ((value << shift) & mix)

  float32 scale = (high - low) * (1.0f / 4294967296.0f);
  unsigned index = m_index;
  while (length--) {
    unsigned index_9  = (index +  9) & 15;
    unsigned index_13 = (index + 13) & 15;
    unsigned index_15 = (index + 15) & 15;

    unsigned state_index = m_state[index];
    unsigned state_index_9 = m_state[index_9];
    unsigned state_index_13 = m_state[index_13];
    unsigned state_index_15 = m_state[index_15];

    unsigned z1 = MUTATE_LEFT(state_index, 16);
    z1 ^= MUTATE_LEFT(state_index_13, 15);

    unsigned z2 = MUTATE_RIGHT(state_index_9, 11);

    unsigned result0 = z1 ^ z2;
    m_state[index] = result0;

    unsigned result1 = MUTATE_LEFT(state_index_15, 2);
    result1 ^= MUTATE_LEFT(z1, 18);;
    result1 ^= z2 << 28;

    result1 ^= MUTATE_LEFT_MIX(result0, 5, 0xda442d24U);

    index = index_15;
    m_state[index] = result1;

    *dst++ = static_cast<float32>(result1) * scale + low;
  }
  m_index = index;

  #undef MUTATE_LEFT
  #undef MUTATE_RIGHT
  #undef MUTATE_LEFT_MIX
}

}  // namespace arraymath
