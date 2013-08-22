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

//------------------------------------------------------------------------------
// This is an implementation of the WELL512 pseudo random number generator. Its
// main characteristics are that it's very fast and has a very long period
// (2^512 - 1 samples).
//
// This is an SSE2 optimized version that in fact implements four parallel
// WELL512 generators that are interleaved to form a single random number
// stream.
//
// This particular implementation is based on an implementation by
// Peter Pettersson: https://github.com/ppettersson/random_well512a_simd/
//------------------------------------------------------------------------------

#include "x86/RandomSSE2.h"

#if defined(AM_USE_X86) && defined(AM_HAS_SSE2)

namespace arraymath {

AM_INLINE
void RandomSSE2::generate4() {
  #define MUTATE_LEFT(value, shift) _mm_xor_si128(value, _mm_slli_epi32(value, shift))
  #define MUTATE_RIGHT(value, shift) _mm_xor_si128(value, _mm_srli_epi32(value, shift))
  #define MUTATE_LEFT_MIX(value, shift, mix) _mm_xor_si128(value, _mm_and_si128(_mm_slli_epi32(value, shift), mix))

  unsigned index_15 = (m_index + 15) & 15;
  __m128i state_index = m_state[m_index];
  __m128i state_index_9 = m_state[(m_index +  9) & 15];
  __m128i state_index_13 = m_state[(m_index + 13) & 15];
  __m128i state_index_15 = m_state[m_index];
  const __m128i kMix = _mm_set1_epi32(0xda442d24);

  __m128i z1 = _mm_xor_si128(MUTATE_LEFT(state_index, 16), MUTATE_LEFT(state_index_13, 15));
  __m128i z2 = MUTATE_RIGHT(state_index_9, 11);
  __m128i result0 = _mm_xor_si128(z1, z2);
  m_state[m_index] = result0;

  __m128i result1 = MUTATE_LEFT(state_index_15, 2);
  result1 = _mm_xor_si128(result1, MUTATE_LEFT(z1, 18));
  result1 = _mm_xor_si128(result1, _mm_slli_epi32(z2, 28));
  result1 = _mm_xor_si128(result1, MUTATE_LEFT_MIX(result0, 5, kMix));
  m_index = index_15;
  m_state[m_index] = result1;

  _mm_store_si128(&m_generated, result1);

  #undef MUTATE_LEFT
  #undef MUTATE_RIGHT
  #undef MUTATE_LEFT_MIX
}

RandomSSE2::RandomSSE2() {
  // Seed the state (64 32-bit integers).
  uint32* state = reinterpret_cast<uint32*>(&m_state[0]);
  uint32 x = *state++ = 5489;
  for (unsigned i = 1; i < 64; ++i) {
    x = *state++ = 1812433253 * (x ^ (x >> 30)) + i;
  }
  m_index = 0;

  // Generate initial batch.
  generate4();
  m_generated_idx = 0;
}

void RandomSSE2::random(float32 *dst, float32 low, float32 high, size_t length) {
  float32 scale = (high - low) * (1.0f / 4294967296.0f);
  uint32* generated = reinterpret_cast<uint32*>(&m_generated);

  // 1) Line up with the batch generator.
  while (m_generated_idx < 4 && length--) {
    *dst++ = static_cast<float32>(generated[m_generated_idx++]) * scale + low;
  }

  // 2) Unrolled core loop.
  while (length >= 4) {
    generate4();
    // TODO(m): Can we SIMDify this somehow?
    *dst++ = static_cast<float32>(generated[0]) * scale + low;
    *dst++ = static_cast<float32>(generated[1]) * scale + low;
    *dst++ = static_cast<float32>(generated[2]) * scale + low;
    *dst++ = static_cast<float32>(generated[3]) * scale + low;
    length -= 4;
  }

  // 3) At this point we're aligned, and need to generate a new batch.
  generate4();
  m_generated_idx = 0;

  // 4) Tail loop (max 3 elements).
  while (length--) {
    *dst++ = static_cast<float32>(generated[m_generated_idx++]) * scale + low;
  }
}

} // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE2
