#include "raytracer.h"
#include "bmp_io.h"
#include "random.h"
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <reducer_opadd.h>

Raytracer::Raytracer(bool soft_shadows, bool direct_illumination,
                     bool global_illumination, bool caustics)
      : _lightSource(NULL)
{
   _root = new SceneDagNode();

   this->soft_shadows = soft_shadows;
   this->direct_illumination = direct_illumination;
   this->global_illumination = global_illumination;
   this->caustics = caustics;
}

Raytracer::~Raytracer()
{
   delete _root;
   delete _lightSource;
}

SceneDagNode* Raytracer::addObject(SceneDagNode* parent, SceneObject* obj)
{
   SceneDagNode* node = new SceneDagNode(obj);
   node->parent = parent;
   node->next = NULL;
   node->child = NULL;

   // Add the object to the parent's child list, this means
   // whatever transformation applied to the parent will also
   // be applied to the child.
   if (parent->child == NULL) {
      parent->child = node;
   } else {
      parent = parent->child;
      while (parent->next != NULL) {
         parent = parent->next;
      }
      parent->next = node;
   }

   return node;
}

LightListNode* Raytracer::addLightSource(LightSource* light)
{
   LightListNode* tmp = _lightSource;
   _lightSource = new LightListNode(light, tmp);
   return _lightSource;
}

void Raytracer::rotate(SceneDagNode* node, char axis, float angle)
{
   Matrix4x4 rotation;
   float toRadian = 2 * M_PI / 360.0;
   int i;

   initMatrix(rotation);
   for (i = 0; i < 2; i++) {
      switch (axis) {
      case 'x':
         rotation[0][0] = 1;
         rotation[1][1] = cos(angle * toRadian);
         rotation[1][2] = -sin(angle * toRadian);
         rotation[2][1] = sin(angle * toRadian);
         rotation[2][2] = cos(angle * toRadian);
         rotation[3][3] = 1;
         break;
      case 'y':
         rotation[0][0] = cos(angle * toRadian);
         rotation[0][2] = sin(angle * toRadian);
         rotation[1][1] = 1;
         rotation[2][0] = -sin(angle * toRadian);
         rotation[2][2] = cos(angle * toRadian);
         rotation[3][3] = 1;
         break;
      case 'z':
         rotation[0][0] = cos(angle * toRadian);
         rotation[0][1] = -sin(angle * toRadian);
         rotation[1][0] = sin(angle * toRadian);
         rotation[1][1] = cos(angle * toRadian);
         rotation[2][2] = 1;
         rotation[3][3] = 1;
         break;
      }
      if (i == 0) {
         mulMatrix(node->trans, node->trans, rotation);
         angle = -angle;
      } else {
         mulMatrix(node->invtrans, rotation, node->invtrans);
      }
   }
}

void Raytracer::translate(SceneDagNode* node, Vector3D trans)
{
   Matrix4x4 translation;

   initMatrix(translation);
   translation[0][3] = trans.v[0];
   translation[1][3] = trans.v[1];
   translation[2][3] = trans.v[2];
   mulMatrix(node->trans, node->trans, translation);
   translation[0][3] = -trans.v[0];
   translation[1][3] = -trans.v[1];
   translation[2][3] = -trans.v[2];
   mulMatrix(node->invtrans, translation, node->invtrans);
   if (node->obj)
      node->obj->ident = false;
}

void Raytracer::scale(SceneDagNode* node, Point3D origin, float factor[3])
{
	Matrix4x4 scale;

	initMatrix(scale);
	scale[0][0] = factor[0];
	scale[0][3] = origin.v[0] - factor[0] * origin.v[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin.v[1] - factor[1] * origin.v[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin.v[2] - factor[2] * origin.v[2];
	mulMatrix(node->trans, node->trans, scale);
	scale[0][0] = 1 / factor[0];
	scale[0][3] = origin.v[0] - 1 / factor[0] * origin.v[0];
	scale[1][1] = 1 / factor[1];
	scale[1][3] = origin.v[1] - 1 / factor[1] * origin.v[1];
	scale[2][2] = 1 / factor[2];
	scale[2][3] = origin.v[2] - 1 / factor[2] * origin.v[2];
	mulMatrix(node->invtrans, scale, node->invtrans);
	if (node->obj)
	{
		node->obj->ident = false;
	}
}

void Raytracer::initInvViewMatrix(Matrix4x4 mat, Point3D eye, Vector3D view,
                                  Vector3D up)
{
   Vector3D w;
   view.normalize();
   up = up - up.dot(view) * view;
   up.normalize();
   w = view.cross(up);

   initMatrix(mat);
   mat[0][0] = w.v[0];
   mat[1][0] = w.v[1];
   mat[2][0] = w.v[2];
   mat[0][1] = up.v[0];
   mat[1][1] = up.v[1];
   mat[2][1] = up.v[2];
   mat[0][2] = -view.v[0];
   mat[1][2] = -view.v[1];
   mat[2][2] = -view.v[2];
   mat[0][3] = eye.v[0];
   mat[1][3] = eye.v[1];
   mat[2][3] = eye.v[2];
}

void Raytracer::flattenEntireScene()
{
   // Flatten scene graph before rendering
   printf("Flattening scene graph... ");
   Matrix4x4 mw, wm;
   initMatrix(mw);
   initMatrix(wm);
   flattenScene(_root, mw, wm);
   printf("Done.\n");
}

void Raytracer::flattenScene(SceneDagNode* node, Matrix4x4 parentmw, Matrix4x4 parentwm)
{
	if(node->obj)
	{
		mulMatrix(node->obj->modelToWorld, parentmw, node->trans);
		mulMatrix(node->obj->worldToModel, node->invtrans, parentwm);

		parentmw = node->obj->modelToWorld;
		parentwm = node->obj->worldToModel;
	}

	SceneDagNode *childPtr = node->child;
	while(childPtr != NULL) {
		flattenScene(childPtr, parentmw, parentwm);
		childPtr = childPtr->next;
	}
}

void Raytracer::traverseEntireScene(Ray3D& ray, bool casting)
{
   traverseScene(_root, ray , casting);
}

void Raytracer::traverseScene(SceneDagNode* node, Ray3D& ray, bool casting)
{
   SceneDagNode *childPtr;

   if (node->obj) {
      // Perform intersection.
      if (casting) {
         if (!node->obj->light)
            node->obj->intersect(ray);
      } else
         node->obj->intersect(ray);
   }
   // Traverse the children.
   childPtr = node->child;
   while (childPtr != NULL) {
      traverseScene(childPtr, ray, casting);
      childPtr = childPtr->next;
   }
}

/* Casts a ray to an object, then to lights, then performs reflect and refract calculations */
void Raytracer::computeShading(Ray3D& ray, int depth, bool getDirectly)
{
   if (!depth) {
      ray.col = default_col;
      return;
   }

   ray.intersection.normal.normalize();

   LightListNode* curLight = _lightSource;
   for (;;) {
      if (curLight == NULL) break;
      // Each light source provides its own shading function.
      curLight->light->shade(ray, getDirectly);
      curLight = curLight->next;
   }

   if (ray.intersection.mat->isDiffuse) {
      ray.col.clamp();
      return;
   }

   float n;

   if (ray.intersection.mat->refractive.max() > 0.01) {

      if (ray.dir.dot(ray.intersection.normal) < 0) {
         n = 1.0 / ray.intersection.mat->refr_index;
      } else {
         ray.intersection.normal = -ray.intersection.normal;
         n = ray.intersection.mat->refr_index;
      }

      float cosI = ray.intersection.normal.dot(ray.dir);
      float sinT2 = n * n * (1.0 - cosI * cosI);

      if (sinT2 < 1.0) {
         Vector3D T = n * ray.dir - (n * cosI + sqrt(1.0 - sinT2)) * ray.intersection.normal;

         Ray3D new_ray(ray.intersection.point, T);
         traverseEntireScene(new_ray, false);

         if (!new_ray.intersection.none)
            computeShading(new_ray, depth - 1, getDirectly);
         else
            new_ray.col = default_col;

         ray.col = ray.col + ray.intersection.mat->refractive * new_ray.col;
      } else {
         // Total internal reflection
         Vector3D R = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
                      ray.intersection.normal;

         Ray3D new_ray(ray.intersection.point, R);
         traverseScene(_root, new_ray, false);

         if (!new_ray.intersection.none)
            computeShading(new_ray, depth - 1, getDirectly);
         else
            new_ray.col = default_col;

         ray.col = ray.col + ray.intersection.mat->refractive * new_ray.col;
      }
   }

   // Raytrace reflections
   if (ray.dir.dot(ray.intersection.normal) < 0) {
      Vector3D R = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
                   ray.intersection.normal;

      Ray3D new_ray(ray.intersection.point, R);
      traverseScene(_root, new_ray, false);

      if (!new_ray.intersection.none)
         computeShading(new_ray, depth - 1, getDirectly);
      else
         new_ray.col = default_col;

      ray.col = ray.col + ray.intersection.mat->specular * new_ray.col;
   }

   ray.col.clamp();
}


void Raytracer::initPixelBuffer()
{
   int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
   _rbuffer = new unsigned char[numbytes];
   _gbuffer = new unsigned char[numbytes];
   _bbuffer = new unsigned char[numbytes];
   for (int i = 0; i < _scrHeight; i++) {
      for (int j = 0; j < _scrWidth; j++) {
         _rbuffer[i*_scrWidth+j] = int(default_col.r * 255);
         _gbuffer[i*_scrWidth+j] = int(default_col.r * 255);
         _bbuffer[i*_scrWidth+j] = int(default_col.r * 255);
      }
   }
}

void Raytracer::flushPixelBuffer(const char *file_name)
{
   bmp_write(file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer);
   delete [] _rbuffer;
   delete [] _gbuffer;
   delete [] _bbuffer;
}

// Sets up ray origin and direction in view space,
// image plane is at z = -1.
Point3D origin(0, 0, 150);

void Raytracer::render_aux(int offW, int offH, int side) 
{
	if(side <= 50)
	{
		int minibound = min(offW+side, _scrWidth);
		int minjbound = min(offH+side, _scrHeight);
		cilk_for(int i = offW; i < minibound; ++i)
		{
			for(int j = offH; j < minjbound; ++j)
			{
				Point3D imagePlane;
				imagePlane.v[2] = 49.9;
				imagePlane.v[1] = (99 * j / float(_scrHeight) - 49.5);
				imagePlane.v[0] = (99 * i / float(_scrWidth) - 49.5);

				// Initialize ray with the proper origin and direction.
				Vector3D dir = imagePlane - origin;
				dir.normalize();

				Ray3D ray(imagePlane, dir);

				traverseScene(_root, ray, false);

				// Don't bother shading if the ray didn't hit anything.
				if (!ray.intersection.none) {
					computeShading(ray, 6, false);
					_rbuffer[j*_scrWidth+i] = int(ray.col.r * 255);
					_gbuffer[j*_scrWidth+i] = int(ray.col.g * 255);
					_bbuffer[j*_scrWidth+i] = int(ray.col.b * 255);
				}
			}
		}
		return;
	}
   cilk_spawn render_aux(offW, offH, side/2);
   cilk_spawn render_aux(offW, offH+side/2, side/2);
   cilk_spawn render_aux(offW+side/2, offH, side/2);
   render_aux(offW+side/2, offH+side/2, side/2);
   cilk_sync;
	return;
}

void Raytracer::render(int width, int height, Point3D eye, Vector3D view,
                       Vector3D up, float fov, const char* fileName)
{
   Matrix4x4 viewToWorld;
   _scrWidth = width;
   _scrHeight = height;

   default_col = Colour(0.4, 0.4, 0.4);

   initPixelBuffer();
   initInvViewMatrix(viewToWorld, eye, view, up);
   clock_t start, end;
   start = clock();

   // Construct a ray for each pixel.
   /* a = i*scrWidth + j */
#define NUMCORES 48
   cilk_spawn render_aux(0, 0, width/2);
   cilk_spawn render_aux(0, width/2, width/2);
   cilk_spawn render_aux(width/2, 0, width/2);
   render_aux(width/2, width/2, width/2);
   cilk_sync;
	end = clock();
	printf("time = %f\n", (end-start)/(float)CLOCKS_PER_SEC);
   flushPixelBuffer(fileName);
}

void Raytracer::printProgress(int percent)
{
   printf("\b\b\b\b");
   cout << percent << "%";
   if (percent < 10) {
      cout << " ";
   }
   if (percent < 100) {
      cout << " ";
   }
   fflush(stdout);
}
