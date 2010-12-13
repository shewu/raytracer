#include "random.h"
#include <cilk_mutex.h>

cilk::mutex Random::lock;

void Random::RandomInit(uint32 seed)
{
   // re-seed generator
   mt[0] = seed;
   for (mti = 1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}


void Random::RandomInitByArray(uint32 seeds[], int length)
{
   // seed by more than 32 bits
   int i, j, k;
   RandomInit(19650218UL);
   if (length <= 0) return;
   i = 1;
   j = 0;
   k = (MERS_N > length ? MERS_N : length);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
      i++;
      j++;
      if (i >= MERS_N) {
         mt[0] = mt[MERS_N-1];
         i = 1;
      }
      if (j >= length) j = 0;
   }
   for (k = MERS_N - 1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {
         mt[0] = mt[MERS_N-1];
         i = 1;
      }
   }
   mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
}

uint32 Random::BRandom()
{
   // generate 32 random bits
   uint32 y;

	 lock.lock();

   if (mti >= MERS_N) {
      // generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1; // lower MERS_R bits
      const uint32 UPPER_MASK = -1L  << MERS_R;      // upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {
                                        0, MERS_A
                                     };

      int kk;
      for (kk = 0; kk < MERS_N - MERS_M; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];
      }

      for (; kk < MERS_N - 1; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
      }

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }

   y = mt[mti++];

	 lock.unlock();

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;
   return y;
}


float Random::Random1()
{
   union {float f;
      uint32 i;
   } convert;
   convert.i = (0x7F << 23) | (((1 << 23) - 1) & BRandom());
   return convert.f - 1.0;
}

float Random::Random2()
{
   union {float f;
      uint32 i;
   } convert;
   convert.i = (0x80 << 23) | (((1 << 23) - 1) & BRandom());
   return convert.f - 3.0;
}

int Random::IRandom(int min, int max)
{
   // output random integer in the interval min <= x <= max
   int r;
   r = int((max - min + 1) * BRandom()) + min; // multiply interval with random and truncate
   if (r > max) r = max;
   if (max < min) return 0x80000000;
   return r;
}

