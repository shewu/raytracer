#ifndef ICACHE_H
#define ICACHE_H

#include <cmath>
#include "util.h"

class ICache
{
	private:
		struct Sample
		{
			Point3D pos;
			Vector3D norm;
			float invR0;
			float r0;
			float tolerance;
			Colour irr;
			bool has_irr;
			Sample *next;

			Sample(Point3D& pos, Vector3D& norm);
			Sample(Point3D& pos, Vector3D& norm, float r0, Colour& irr, float tolerance);
			~Sample();
			float weight(Sample& x);
			Colour getIrradiance(Sample& x);
		};

		Sample *first;
		float tolerance;
		float invTolerance;
		float minSpacing;
		float maxSpacing;
		float find(Sample& x);

		friend struct ICache::Sample;

	public:

		ICache(float tolerance, float minSpacing);
		~ICache();
		void insert(Point3D& p, Vector3D& n, float r0, Colour& irr);
		bool getIrradiance(Point3D& p, Vector3D& n, Colour *c);
};

#endif
