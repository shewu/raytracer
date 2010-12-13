#include <cmath>
#include <iostream>

#include "config.h"
#include "scene_object.h"

#define UPDATE_RAY(N1, N2, N3, MAT)\
    ray.intersection.point = p.transform(modelToWorld);\
    ray.intersection.t_value = lambda;\
    ray.intersection.normal = Vector3D(N1, N2, N3).transform(modelToWorld);\
    ray.intersection.none = false;\
    ray.intersection.mat = MAT

#define UPDATE_RAY_IDENT(N1, N2, N3, MAT, lambda, p)\
    ray.intersection.point = p;\
    ray.intersection.t_value = lambda;\
    ray.intersection.normal = Vector3D(N1, N2, N3);\
    ray.intersection.none = false;\
    ray.intersection.mat = MAT

bool Square::intersect(Ray3D& ray)
{

   Point3D origin = ray.origin.transform(worldToModel);
   Vector3D dir = ray.dir.transform(worldToModel);

   float lambda = (49.99 - origin.v[1]) / dir.v[1];
   if (lambda > EPSILON) {
      Point3D p = origin + lambda * dir;
      if (p.v[0] <= 16 && p.v[0] >= -16 && p.v[2] <= 16 && p.v[2] >= -16) {
         UPDATE_RAY(0, -1, 0, mat);

         return true;
      }
   }

   return false;
}

bool Cube::intersect(Ray3D& ray)
{
   Point3D p;
   float lambda;
   Point3D origin = ray.origin;
   Vector3D dir = ray.dir;

   lambda = (50 - origin.v[2]) / dir.v[2];
   if (lambda > EPSILON && lambda < ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[0] <= 50 && p.v[0] >= -50 && p.v[1] <= 50 && p.v[1] >= -50) {
         UPDATE_RAY_IDENT(0, 0, d, mat1, lambda, p);

         return true;
      }
   }

   lambda = (-50 - origin.v[2]) / dir.v[2];
   if (lambda > EPSILON && lambda < ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[0] <= 50 && p.v[0] >= -50 && p.v[1] <= 50 && p.v[1] >= -50) {
         UPDATE_RAY_IDENT(0, 0, -d, mat2, lambda, p);

         return true;
      }
   }

   lambda = (50 - origin.v[1]) / dir.v[1];
   if (lambda > EPSILON && lambda < ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[0] <= 50 && p.v[0] >= -50 && p.v[2] <= 50 && p.v[2] >= -50) {
         UPDATE_RAY_IDENT(0, d, 0, mat3, lambda, p);

         return true;
      }
   }

   lambda = (-50 - origin.v[1]) / dir.v[1];
   if (lambda > EPSILON && lambda < ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[0] <= 50 && p.v[0] >= -50 && p.v[2] <= 50 && p.v[2] >= -50) {
         UPDATE_RAY_IDENT(0, -d, 0, mat4, lambda, p);

         return true;
      }
   }

   lambda = (50 - origin.v[0]) / dir.v[0];
   if (lambda > EPSILON && lambda <= ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[1] <= 50 && p.v[1] >= -50 && p.v[2] <= 50 && p.v[2] >= -50) {
         UPDATE_RAY_IDENT(d, 0, 0, mat5, lambda, p);

         return true;
      }
   }

   lambda = (-50 - origin.v[0]) / dir.v[0];
   if (lambda > EPSILON && lambda <= ray.intersection.t_value) {
      p = origin + lambda * dir;
      if (p.v[1] <= 50 && p.v[1] >= -50 && p.v[2] <= 50 && p.v[2] >= -50) {
         UPDATE_RAY_IDENT(-d, 0, 0, mat6, lambda, p);

         return true;
      }
   }

   return false;
}

bool Sphere::intersect(Ray3D& ray)
{
   Point3D p;
   float lambda;
   float A, B, C, D;
   Point3D origin = ray.origin.transform(worldToModel);
   Vector3D dir = ray.dir.transform(worldToModel);

   float radius_sq = 20 * 20;
   A = dir.dot(dir);
   B = 2 * origin.dot(dir);
   C = origin.dot(origin) - radius_sq;
   D = B * B - (4 * A * C);

   if (D < 0) return false;
   else if (D > 0) {
      lambda = (-B - sqrt(D)) / (2 * A);

      if (lambda > EPSILON && lambda < ray.intersection.t_value) {
         p = origin + lambda * dir;
         UPDATE_RAY(p.v[0], p.v[1], p.v[2], mat);

         return true;
      }

      lambda = (-B + sqrt(D)) / (2 * A);

      if (lambda > EPSILON && lambda < ray.intersection.t_value) {
         p = origin + lambda * dir;
         UPDATE_RAY(p.v[0], p.v[1], p.v[2], mat);

         return true;
      }
   } else {
      lambda = -B / (2 * A);

      if (lambda > EPSILON && lambda < ray.intersection.t_value) {
         p = origin + lambda * dir;
         UPDATE_RAY(p.v[0], p.v[1], p.v[2], mat);

         return true;
      }
   }

   return false;
}


bool DisplacedSurface::checkIntersectionGrid(
   int xcoord, int zcoord,
   Point3D &origin, Vector3D &dir,
   Ray3D &ray,
   Matrix4x4 modelToWorld)
{

   // Each square grid is made up of to triages (A and B).  This code checks
   // to see if the ray intersections with either of them.  If the ray
   // intersects with both triangles, we determine which is closer.

   /*     a         d
    *       + -----+
    *       |     /|
    *       | A  / |
    *       |   /  |
    *       |  /   |
    *       | /  B |
    *       |/     |
    *       +------+
    *     b         c
    */

   int i = -1;
   Point3D x1, x2;
   Vector3D n1, n2;
   float t1, t2;
   Point3D a, b, c, d;

   a = points[xcoord + 1][zcoord];
   b = points[xcoord][zcoord];
   c = points[xcoord][zcoord + 1];
   d = points[xcoord + 1][zcoord + 1];

   n1 = (b - a).cross(d - a);

   t1 = (a - origin).dot(n1) / (dir.dot(n1));
   x1 = origin + t1 * dir;
   if (((b - a).cross(x1 - a)).dot(n1) >= 0 &&
       ((d - b).cross(x1 - b)).dot(n1) >= 0 &&
       ((a - d).cross(x1 - d)).dot(n1) >= 0) {
      i = 1;
   }

   n2 = (d - c).cross(b - c);

   t2 = (c - origin).dot(n2) / (dir.dot(n2));
   if (i != 1 || t2 < t1) {
      x2 = origin + t2 * dir;
      if (((d - c).cross(x2 - c)).dot(n2) >= 0 &&
          ((b - d).cross(x2 - d)).dot(n2) >= 0 &&
          ((c - b).cross(x2 - b)).dot(n2) >= 0)
         i = 2;
   }

   if (i == 1 && t1 > EPSILON && t1 < ray.intersection.t_value) {
      n1.normalize();
      UPDATE_RAY_IDENT(n1.v[0], n1.v[1], n1.v[2], mat, t1, x1);
      return true;
   }
   if (i == 2 && t2 > EPSILON && t2 < ray.intersection.t_value) {
      n2.normalize();
      UPDATE_RAY_IDENT(n2.v[0], n2.v[1], n2.v[2], mat, t2, x2);
      return true;
   }

   return false;
}

bool DisplacedSurface::intersect(Ray3D& ray)
{
   Point3D origin = ray.origin;
   Vector3D dir = ray.dir;

   bool intersected = false;

   // Check each triangle in grid to see if it intersects with this ray
   for (int x1 = 0; x1 < xcoords; x1++) {
      for (int z1 = 0; z1 < zcoords; z1++) {
         intersected |= checkIntersectionGrid(x1, z1, origin, dir, ray,
                                              modelToWorld);
      }
   }

   return intersected;
}

void DisplacedSurface::makeSurface()
{
   /* Extents are -50 to +50 in x-z plane, offset by height */
   int z = 0;
   int x = 0;
   float time = 0;

   for (x = 0; x < xcoords + 1; x++) {
      for (z = 0; z < zcoords + 1; z++) {

         float xc = (x * 100.0) / xcoords - 50.0;
         float zc = (z * 100.0) / zcoords - 50.0;
         float y = perlin_noise.Get(xc / 1.2, zc, time) + middle;

         points[x][z] = Point3D(xc, y, zc);
      }
   }

}
