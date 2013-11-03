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
// This is a NEON optimized version that in fact implements four parallel
// WELL512 generators that are interleaved to form a single random number
// stream.
//
// This particular implementation is based on an implementation by
// Peter Pettersson: https://github.com/ppettersson/random_well512a_simd/
//------------------------------------------------------------------------------

#include "arm/RandomNEON.h"

#if defined(AM_USE_ARM) && defined(AM_HAS_NEON)

#include <arm_neon.h>

#include <algorithm>

#include "Random.h"

namespace arraymath {

class RandomNEON : public Random {
 public:
  virtual ~RandomNEON();

  virtual void random(float32 *dst, float32 low, float32 high, size_t length);

  bool init();

 private:
  float32x4_t generate4();

  uint32x4_t* m_state;
  float32* m_generated;

  unsigned m_index;
  unsigned m_generated_idx;
};

#define mutateLeft(value, shift) \
  veorq_u32(value, vshlq_n_u32(value, shift))

#define mutateRight(value, shift) \
  veorq_u32(value, vshrq_n_u32(value, shift))

#define mutateLeftMix(value, shift, mix) \
  veorq_u32(value, vandq_u32(vshlq_n_u32(value, shift), mix))

AM_INLINE
float32x4_t RandomNEON::generate4() {
  unsigned index_15 = (m_index + 15) & 15;
  uint32x4_t state_index = m_state[m_index];
  uint32x4_t state_index_9 = m_state[(m_index +  9) & 15];
  uint32x4_t state_index_13 = m_state[(m_index + 13) & 15];
  uint32x4_t state_index_15 = m_state[m_index];
  const uint32x4_t kMix = vdupq_n_u32(0xda442d24);

  uint32x4_t z1 = veorq_u32(mutateLeft(state_index, 16), mutateLeft(state_index_13, 15));
  uint32x4_t z2 = mutateRight(state_index_9, 11);
  uint32x4_t result0 = veorq_u32(z1, z2);
  m_state[m_index] = result0;

  uint32x4_t result1 = mutateLeft(state_index_15, 2);
  result1 = veorq_u32(result1, mutateLeft(z1, 18));
  result1 = veorq_u32(result1, vshlq_n_u32(z2, 28));
  result1 = veorq_u32(result1, mutateLeftMix(result0, 5, kMix));
  m_index = index_15;
  m_state[m_index] = result1;

  // Convert the integers to floating point numbers in the range [0.5, 1.0)
  // using some bit trickery.
  const uint32x4_t kExponentMask = vdupq_n_u32(0x007fffff);
  const uint32x4_t kExponent = vdupq_n_u32(0x3f000000);
  result1 = vorrq_u32(vandq_u32(result1, kExponentMask), kExponent);
  return vreinterpretq_f32_u32(result1);
}

bool RandomNEON::init() {
  // Allocate memory for the NEON state and the 4-float batch buffer.
  // Note: They have to be 16-byte aligned, which they might not be if we put
  // it directly in the class and we're compiling for a 32-bit architecture,
  // which only guarantees 8-byte struct alignment.
  m_state = new uint32x4_t[16 + 1];
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
  vst1q_f32(m_generated, generate4());
  m_generated_idx = 0;

  return true;
}

RandomNEON::~RandomNEON() {
  delete[] m_state;
}

void RandomNEON::random(float32 *dst, float32 low, float32 high, size_t length) {
  // Scaling factors for turning the range [0.5, 1.0) to [low, high).
  float32 scale = (high - low) * 2.0f;
  float32 offset = low - 0.5f * scale;

  // 1) Line up with the batch generator.
  unsigned numLineUp = std::min(m_generated_idx ? 4 - m_generated_idx : 0, static_cast<unsigned>(length));
  length -= numLineUp;
  while (numLineUp--) {
    *dst++ = m_generated[m_generated_idx++] * scale + offset;
  }
  if (AM_UNLIKELY(m_generated_idx < 4)) {
    // If we didn't finish the current batch, exit now (we'll continue the line
    // up operation during the next call).
    return;
  }

  // 2) Unrolled core loop.
  float32x4_t _scale = vdupq_n_f32(scale);
  float32x4_t _offset = vdupq_n_f32(offset);
  while (length >= 4) {
    float32x4_t rnd = generate4();
    vst1q_f32(dst, vaddq_f32(vmulq_f32(rnd, _scale), _offset));
    dst += 4;
    length -= 4;
  }

  // 3) At this point we need to generate a new batch.
  vst1q_f32(m_generated, generate4());
  m_generated_idx = 0;

  // 4) Tail loop (max 3 elements).
  while (length--) {
    *dst++ = m_generated[m_generated_idx++] * scale + offset;
  }
}

Random* RandomNEONFactory::create() {
  RandomNEON* random = new RandomNEON();
  if (random && !random->init()) {
    delete random;
    return NULL;
  }
  return random;
}

} // namespace arraymath

#endif // AM_USE_ARM && AM_HAS_NEON
