#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <cilk_mutex.h>

typedef   signed long int32;
typedef unsigned long uint32;

class Random                  // encapsulate random number generator
{
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
public:
	 Random(uint32 seed=2)
	 {
			RandomInit(seed);
	 }
   void RandomInit(uint32 seed);        // re-seed
   void RandomInitByArray(uint32 seeds[], int length); // seed by more than 32 bits
   int IRandom(int min, int max);       // output random integer
   float Random1();                     // output random float
   float Random2();                     // output random float
   uint32 BRandom();                    // output random bits
	 static cilk::mutex lock;
private:
   uint32 mt[MERS_N];                   // state vector
   int mti;                             // index into mt
};

/*
class Random						// HyperObject for a random number generator
{
	struct RandMonoid: cilk::monoid_base<RandomBase>
	{
		void reduce(RandomBase *L, RandomBase *R) {}
		void identity(RandomBase* r) const
		{
			uint32 i;
			new (r) RandomBase(i);
		}
	};
	cilk::reducer<RandMonoid> reducer;

public:
	void RandomInit(uint32 seed) { reducer().RandomInit(seed); }
	int IRandom(int min, int max) { return reducer().IRandom(min, max); }
	float Random1() { return reducer().Random1(); }
	float Random2() { return reducer().Random2(); }
	uint32 BRandom() { return reducer().BRandom(); }
};*/

#endif
