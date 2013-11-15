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

#include <emmintrin.h>

#include <algorithm>

#include "Random.h"

namespace arraymath {

class RandomSSE2 : public Random {
 public:
  virtual ~RandomSSE2();

  virtual void random(float32 *dst, float32 low, float32 high, size_t length);

  bool init();

 private:
  __m128 generate4();

  __m128i* m_state;
  float32* m_generated;

  unsigned m_index;
  unsigned m_generated_idx;
};

namespace {

__m128i mutateLeft(__m128i value, int shift) {
  return _mm_xor_si128(value, _mm_slli_epi32(value, shift));
}

__m128i mutateRight(__m128i value, int shift) {
  return _mm_xor_si128(value, _mm_srli_epi32(value, shift));
}

__m128i mutateLeftMix(__m128i value, int shift, __m128i mix) {
  return _mm_xor_si128(value, _mm_and_si128(_mm_slli_epi32(value, shift), mix));
}

} // anonymous namespace

AM_INLINE
__m128 RandomSSE2::generate4() {
  unsigned index_15 = (m_index + 15) & 15;
  __m128i state_index = m_state[m_index];
  __m128i state_index_9 = m_state[(m_index +  9) & 15];
  __m128i state_index_13 = m_state[(m_index + 13) & 15];
  __m128i state_index_15 = m_state[m_index];
  const __m128i kMix = _mm_set1_epi32(0xda442d24);

  __m128i z1 = _mm_xor_si128(mutateLeft(state_index, 16), mutateLeft(state_index_13, 15));
  __m128i z2 = mutateRight(state_index_9, 11);
  __m128i result0 = _mm_xor_si128(z1, z2);
  m_state[m_index] = result0;

  __m128i result1 = mutateLeft(state_index_15, 2);
  result1 = _mm_xor_si128(result1, mutateLeft(z1, 18));
  result1 = _mm_xor_si128(result1, _mm_slli_epi32(z2, 28));
  result1 = _mm_xor_si128(result1, mutateLeftMix(result0, 5, kMix));
  m_index = index_15;
  m_state[m_index] = result1;

  // Convert the integers to floating point numbers in the range [0.5, 1.0)
  // using some bit trickery.
  const __m128i kExponentMask = _mm_set1_epi32(0x007fffff);
  const __m128i kExponent = _mm_set1_epi32(0x3f000000);
  result1 = _mm_or_si128(_mm_and_si128(result1, kExponentMask), kExponent);
  return _mm_castsi128_ps(result1);
}

bool RandomSSE2::init() {
  // Allocate memory for the SSE state and the 4-float batch buffer.
  // Note: They have to be 16-byte aligned, which they might not be if we put
  // it directly in the class and we're compiling for a 32-bit architecture,
  // which only guarantees 8-byte struct alignment.
  m_state = new __m128i[16 + 1];
  if (!m_state) {
    return false;
  }
  m_generated = reinterpret_cast<float32*>(&m_state[16]);

  // Seed the state (64 32-bit integers).
  uint32* state = reinterpret_cast<uint32*>(m_state);
  uint32 x = *state++ = 5489;
  for (unsigned i = 1; i < 64; ++i) {
    x = *state++ = 1812433253 * (x ^ (x >> 30)) + i;
  }
  m_index = 0;

  // Generate initial batch.
  _mm_store_ps(m_generated, generate4());
  m_generated_idx = 0;

  return true;
}

RandomSSE2::~RandomSSE2() {
  delete[] m_state;
}

void RandomSSE2::random(float32 *dst, float32 low, float32 high, size_t length) {
  // Scaling factors for turning the range [0.5, 1.0) to [low, high).
  float32 scale = (high - low) * 2.0f;
  float32 offset = low - 0.5f * scale;

  // 1) Line up with the batch generator.
  unsigned numLineUp = std::min(m_generated_idx ? 4 - m_generated_idx : 0, static_cast<unsigned>(length));
  length -= numLineUp;
  while (numLineUp--) {
    *dst++ = m_generated[m_generated_idx++] * scale + offset;
  }
  if (AM_UNLIKELY(m_generated_idx > 0 && m_generated_idx < 4)) {
    // If we didn't finish the current batch, exit now (we'll continue the line
    // up operation during the next call).
    return;
  }

  // 2) Unrolled core loop.
  __m128 _scale = _mm_set1_ps(scale);
  __m128 _offset = _mm_set1_ps(offset);
  while (length >= 4) {
    __m128 rnd = generate4();
    _mm_storeu_ps(dst, _mm_add_ps(_mm_mul_ps(rnd, _scale), _offset));
    dst += 4;
    length -= 4;
  }

  // 3) At this point we need to generate a new batch.
  _mm_store_ps(m_generated, generate4());
  m_generated_idx = 0;

  // 4) Tail loop (max 3 elements).
  while (length--) {
    *dst++ = m_generated[m_generated_idx++] * scale + offset;
  }
}

Random* RandomSSE2Factory::create() {
  RandomSSE2* random = new RandomSSE2();
  if (random && !random->init()) {
    delete random;
    return NULL;
  }
  return random;
}

} // namespace arraymath

#endif // AM_USE_X86 && AM_HAS_SSE2
