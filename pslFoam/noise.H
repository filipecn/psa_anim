// functions extracted from both sources
// https://github.com/AcademySoftwareFoundation/OpenImageIO/blob/master/src/include/OpenImageIO/hash.h
// https://github.com/AcademySoftwareFoundation/OpenShadingLanguage/blob/main/src/include/OSL/oslnoise.h

uint32_t rotl32(uint32_t x, int k)
{
  return (x << k) | (x >> (32 - k));
}

// Mix up the bits of a, b, and c (changing their values in place).
inline void bjmix (uint32_t &a, uint32_t &b, uint32_t &c)
{
  a -= c;  a ^= rotl32(c, 4);  c += b;
  b -= a;  b ^= rotl32(a, 6);  a += c;
  c -= b;  c ^= rotl32(b, 8);  b += a;
  a -= c;  a ^= rotl32(c,16);  c += b;
  b -= a;  b ^= rotl32(a,19);  a += c;
  c -= b;  c ^= rotl32(b, 4);  b += a;
}

// Mix up and combine the bits of a, b, and c (doesn't change them, but
// returns a hash of those three original values).  21 ops
inline uint32_t bjfinal (uint32_t a, uint32_t b, uint32_t c=0xdeadbeef)
{
  c ^= b; c -= rotl32(b,14);
  a ^= c; a -= rotl32(c,11);
  b ^= a; b -= rotl32(a,25);
  c ^= b; c -= rotl32(b,16);
  a ^= c; a -= rotl32(c,4);
  b ^= a; b -= rotl32(a,14);
  c ^= b; c -= rotl32(b,24);
  return c;
}


float bits_to_01 (unsigned int bits) {
  // divide by 2^32-1
  // Calculate inverse constant with double precision to avoid
  //     warning: implicit conversion from 'unsigned int' to 'float' changes value from 4294967295 to 4294967296
  constexpr float convertFactor = static_cast<float>(static_cast<double>(1.0) / static_cast<double>(std::numeric_limits<unsigned int>::max()));
  return bits * convertFactor;
}

unsigned int inthash (const unsigned int k0, const unsigned int k1, const unsigned int k2, const unsigned int k3) {
  // now hash the data!
  unsigned int start_val = 0xdeadbeef + (4 << 2) + 13;

  unsigned int a = start_val + k0;
  unsigned int b = start_val + k1;
  unsigned int c = start_val + k2;
  bjmix(a, b, c);
  a += k3;
  c = bjfinal(a, b, c);
  return c;
}

float cellnoise(vector p, float f) {
  return bits_to_01(inthash(
        floor(p.x()),
        floor(p.y()),
        floor(p.z()),
        floor(f)));
}

float clamp(float x, float a, float b) {
  return min(b, max(x, a));
}

float smoothstep(float a, float b, float x) {
  float t = clamp((x - a) / (b - a), 0., 1.);
  return t * t * (3 - 2 * t);
}

double distance(const vector& a, const vector& b) {
  return mag(a - b);
}

/*
 * Copyright (c) 2016
 *  Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 */
  
/// Alligator Noise is provided by Side Effects Software Inc. and is licensed
/// under a Creative Commons Attribution-ShareAlike 4.0 International License.
///
/// Author:  "Side Effects Software Inc"
/// Source:  "http://www.sidefx.com/docs/hdk15.0/alligator_2alligator_8_c-example.html"
/// License: "http://creativecommons.org/licenses/by-sa/4.0/"
///
/// Translated and modified by Ivan Mavrov, Chaos Group Ltd. 2016
/// Contact: ivan.mavrov@chaosgroup.com

 
/// 3D Alligator noise implementation.
/// Returned values are in the [0, 1] range.
float alligatorNoise3D(vector position) {
  vector cellOffsets[27] = {
    vector( 0,  0,  0),
    vector( 1,  0,  0),
    vector( 1,  1,  0),
    vector( 0,  1,  0),
    vector(-1,  1,  0),
    vector(-1,  0,  0),
    vector(-1, -1,  0),
    vector( 0, -1,  0),
    vector( 1, -1,  0),

    vector( 0,  0, -1),
    vector( 1,  0, -1),
    vector( 1,  1, -1),
    vector( 0,  1, -1),
    vector(-1,  1, -1),
    vector(-1,  0, -1),
    vector(-1, -1, -1),
    vector( 0, -1, -1),
    vector( 1, -1, -1),

    vector( 0,  0,  1),
    vector( 1,  0,  1),
    vector( 1,  1,  1),
    vector( 0,  1,  1),
    vector(-1,  1,  1),
    vector(-1,  0,  1),
    vector(-1, -1,  1),
    vector( 0, -1,  1),
    vector( 1, -1,  1)
  };

  vector iPosition(
      floor(position.x()),
      floor(position.y()),
      floor(position.z()));

  float firstReverseSmoothPointDistance = 0.0;
  float secondReverseSmoothPointDistance = 0.0;

  for (int cellIndex = 0; cellIndex < 27; ++cellIndex) {
    point cellOrigin = iPosition + cellOffsets[cellIndex];
    float cn = cellnoise(cellOrigin, 0.0);
    vector cellPointOffset(cn,cn,cn) ;
    point cellPointPosition = cellOrigin + cellPointOffset;

    float cellPointDistance = distance(position, cellPointPosition);

    if (cellPointDistance < 1.0) {
      float reverseSmoothDistance = smoothstep(0.0, 1.0, 1.0 - cellPointDistance);

      float distanceMultiplier = float(cellnoise(cellOrigin, 1.0));
      reverseSmoothDistance *= distanceMultiplier;

      if (firstReverseSmoothPointDistance < reverseSmoothDistance) {
        secondReverseSmoothPointDistance = firstReverseSmoothPointDistance;
        firstReverseSmoothPointDistance = reverseSmoothDistance;
      } else {
        if (secondReverseSmoothPointDistance < reverseSmoothDistance) {
          secondReverseSmoothPointDistance = reverseSmoothDistance;
        }
      }
    }
  }

  return firstReverseSmoothPointDistance - secondReverseSmoothPointDistance;
}

/// 3D Fractal Alligator noise implementation.
/// Returned values are in the [-1, 1] range.
float fractalAlligatorNoise3D(
    point position,
    float lacunarity, // Houdini 2.0
    float gain,       // Houdini rough
    int octaveCount   // Houdini turbulence - 1
    ) {
  float noiseValue = 0.0;

  float amplitude = 1.0;

  for (int octave = 0; octave < octaveCount; ++octave) {
    noiseValue += amplitude * (alligatorNoise3D(position) - 0.5);      
    position *= lacunarity;
    amplitude *= gain;
  }

  return noiseValue;
}

float fractalAllicatorNoise(
    vector P,
    // "Initial sampling position multiplier that affects the overall granularity."
    float start_frequency = 10.0,
    // "Offsets the initial sampling position effectively shifting the pattern in the specified direction."
    vector start_offset = vector(0, 0, 0),
    // "Position (frequency) multiplier per iteration."
    float lacunarity = 2.0,
    //[[ string description = "Amplitude multiplier per iteration."
    float gain = 0.5,
    //[[ string description = "Number of fractal iterations."
    int octaves = 8,
    // "The power of the falloff applied to the final result."
    float attenuation = 1.0
    ) {
  point startPosition = start_frequency * P - start_offset;

  float noiseValue = fractalAlligatorNoise3D(startPosition, lacunarity, gain, octaves);

  noiseValue = 0.5 * noiseValue + 0.5;
  noiseValue = pow(noiseValue, attenuation);

  return noiseValue;
}
