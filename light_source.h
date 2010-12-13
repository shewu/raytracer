#ifndef LIGHT_SOURCE_H
#define LIGHT_SOURCE_H

#include "util.h"
#include "icache.h"
#include "random.h"
#include "photonmap.h"

class Raytracer;

// Base class for a light source.
class LightSource
{
public:
   virtual void shade(Ray3D&, bool getDirectly) = 0;
   virtual ~LightSource()
   { };
};

// Square photon mapping light
class SquarePhotonLight : public LightSource
{
private:
   BalancedPhotonMap *bmap;
   BalancedPhotonMap *cmap;
   Colour col;
   Random r;
   Raytracer *raytracer;
   ICache icache;
   Colour light_col;
   void globalIllumination(Ray3D& ray, bool getDirectly);
   void causticIllumination(Ray3D& ray);
   void directIllumination(Ray3D& ray);

public:

   SquarePhotonLight(Colour col, Raytracer *raytracer);
   ~SquarePhotonLight();

   void initTransformMatrix(Matrix4x4 mat, Vector3D& w);
   void shade(Ray3D& ray, bool getDirectly);
   Vector3D getRandLambertianDir(Vector3D& normal);
   void tracePhotons(int num, int castics_num);
};

#endif
