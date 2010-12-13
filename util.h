#ifndef UTIL_H
#define UTIL_H

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#ifndef NULL
#define NULL 0
#endif

#include <float.h>
#include <emmintrin.h>
#include <xmmintrin.h>

typedef float Matrix4x4[4][4] __attribute__ ((aligned(16)));

inline float rsqrtss_sqrt(float x)
{
	/* use rsqrtss * x */
	__m128 a = _mm_load_ss(&x);
	float ret[1] = {0};
	_mm_store_ss(ret, _mm_mul_ss(a, _mm_rsqrt_ss(a)));
	return *ret;
}

// Initializes mat to the identity matrix.
void initMatrix(Matrix4x4 mat);

// Implements ret = mat1*mat2 .
void mulMatrix(Matrix4x4 ret, Matrix4x4 mat1, Matrix4x4 mat2);

struct Vector3D
{
   float v[4] __attribute__ ((aligned(16)));

   Vector3D()
   {
	   v[0] = 0;
	   v[1] = 0;
	   v[2] = 0;
	   v[3] = 0;
   }
   Vector3D(float x, float y, float z)
   {
	   v[0] = x;
	   v[1] = y;
	   v[2] = z;
	   v[3] = 0;
   }
   inline float dot(const Vector3D& v) const
   {
	   float result[4] = {0,0,0,0};
	   /* seems like adding memory to the clobber list is faster than
		* declaring a volatile array */
	   __asm __volatile__(".intel_syntax noprefix\n\t"
			   "movaps xmm0, [rbx] \n\t"
			   "movaps xmm1, [rcx] \n\t"
			   "mulps xmm0, xmm1 \n\t"
			   "xorps xmm2, xmm2 \n\t"
			   "haddps xmm2, xmm0 \n\t"
			   "haddps xmm2, xmm2 \n\t"
			   "movaps [rax], xmm2 \n\t"
			   ".att_syntax prefix \n\t"
			   : : "c" (this->v), "b" (v.v), "a" (result) : "xmm0", "xmm1", "xmm2", "memory");
	   /* aha, so the first time we access result[0], it's 0, but the second
		* time we do so, the value appears o_o. however, assignment ops don't
		* solve the problem. */
	   /* how is the function returning the correct value when it can't print
		* the value in the line above?! looks like printf is going nuts */
	   return *(result+1);
   }
   Vector3D cross(const Vector3D& v) const;
   Vector3D operator+(const Vector3D& v) const;
   Vector3D operator-(const Vector3D& v) const;
   Vector3D operator-() const;
   float normalize();
   float mag();

   // Returns a transformed vector by right-multiplying it to mat.
   Vector3D transform(Matrix4x4 mat) const;

   // Returns a transformed vector by right-multiplying it to transpose(mat),
   // this is useful for transforming normals back to world space.
   Vector3D transformAsNormal(Matrix4x4 mat) const;
};

// Scalar multiplication.
Vector3D operator*(float s, const Vector3D& v);

struct Point3D
{
   float v[4] __attribute__ ((aligned(16)));
   Point3D()
   {
	   v[0] = 0;
	   v[1] = 0;
	   v[2] = 0;
	   v[3] = 0;
   }
   Point3D(float x, float y, float z) 
   {
	   v[0] = x;
	   v[1] = y;
	   v[2] = z;
	   v[3] = 0;
   }
   Point3D operator+(const Vector3D& v) const;

   Point3D operator-(const Vector3D& v) const;
   /* holy shit SSEing this speeds up the function by 10x */
   inline Vector3D operator-(const Point3D& p) const
   {
	   Vector3D ret;
	   __m128 thisv = _mm_load_ps(this->v);
	   __m128 pv = _mm_load_ps(p.v);
	   __m128 retv = _mm_sub_ps(thisv, pv);
	   _mm_store_ps(ret.v, retv);
	   return ret;
   }
   float operator[](int ind) const;
   float &operator[](int ind);
   Point3D transform(Matrix4x4 mat) const;
   float dot(const Point3D& v) const;
   float distanceTo(const Point3D& p) const;
   //float distance2To(const Point3D& p) const;
   float dot(const Vector3D& v) const;
};

struct Colour
{
   Colour() : r(0), g(0), b(0)
   {}
   Colour(float c) : r(c), g(c), b(c)
   {}
   Colour(float r, float g, float b) : r(r), g(g), b(b)
   {}
   Colour operator*(const Colour& c) const;
   Colour operator+(const Colour& c) const;
   Colour operator-(const Colour& c) const;
   Colour operator+(const float f) const;
   void operator+=(const Colour& c);
   void operator*=(const Colour& c);
   float max();

   // clamp the colour channels to 1.0
   void clamp();

   float r;
   float g;
   float b;
};

Colour operator*(float s, const Colour& c);
Colour operator/(const Colour& c, float s);

struct Material
{
   Material(Colour diffuse, Colour specular, float exp,
            Colour refractive, float refr_index, bool isDiffuse,
            bool isSpecular = false, bool light = false) :
         diffuse(diffuse), specular(specular), refractive(refractive),
         specular_exp(exp), refr_index(refr_index), isDiffuse(isDiffuse),
         isSpecular(isSpecular), light(light)
   {}
   // Probability of diffusion .
   Colour diffuse;
   // Probability of reflection.
   Colour specular;
   // Probability of refraction.
   Colour refractive;
   // Specular exponent (for Phong shading).
   float specular_exp;
   // Material colour.
   Colour col;
   // Refractive index.
   float refr_index;

   bool isDiffuse;
   bool isSpecular;
   bool light;
};

struct Intersection
{
   Intersection() : none(true), t_value(FLT_MAX)
   {}
   // Location of intersection.
   Point3D point;
   // Normal at the intersection.
   Vector3D normal;
   // Material at the intersection.
   Material* mat;
   // This is used when you need to intersect multiply objects and
   // Set to true when no intersection has occurred.
   bool none;
   // Distance to nearest intersection.
   float t_value;
   // Intersected with light?
   bool light;
};

// Ray structure.
struct Ray3D
{
   Ray3D()
   {}
   Ray3D(Point3D p, Vector3D v) : origin(p), dir(v), col(0, 0, 0)
   {}
   // Origin and direction of the ray.
   Point3D origin;
   Vector3D dir;
   // Intersection status, should be computed by the intersection
   // function.
   Intersection intersection;
   // Current colour of the ray, should be computed by the shading
   // function.
   Colour col;
};


#endif
