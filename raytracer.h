#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "util.h"
#include "scene_object.h"
#include "light_source.h"

#include <vector>
using namespace std;

// Linked list containing light sources in the scene.
struct LightListNode
{
   LightListNode() : light(NULL), next(NULL)
   {}
   LightListNode(LightSource* light, LightListNode* next = NULL) :
         light(light), next(next)
   {}
   ~LightListNode()
   {
      delete light;
      delete next;
   }
   LightSource* light;
   LightListNode* next;
};

// The scene graph, containing objects in the scene.
struct SceneDagNode
{
   SceneDagNode() :
         obj(NULL), next(NULL), parent(NULL), child(NULL)
   {
      initMatrix(trans);
      initMatrix(invtrans);
   }

   SceneDagNode(SceneObject* obj) :
         obj(obj), next(NULL), parent(NULL), child(NULL)
   {
      initMatrix(trans);
      initMatrix(invtrans);
   }

   ~SceneDagNode()
   {
      delete obj;
      delete child;
      delete next;
   }

   // Pointer to geometry primitive, used for intersection.
   SceneObject* obj;
   // Pointer to material of the object, used in shading.
   // Each node maintains a transformation matrix, which maps the
   // geometry from object space to world space and the inverse.
   Matrix4x4 trans;
   Matrix4x4 invtrans;

   // Internal structure of the tree, you shouldn't have to worry
   // about them.
   SceneDagNode* next;
   SceneDagNode* parent;
   SceneDagNode* child;
};

class Raytracer
{
public:
   Raytracer(bool soft_shadows = true, bool direct_illumination = true,
             bool global_illumination = true, bool caustics = true);
   ~Raytracer();

   // Renders an image fileName with width and height and a camera
   // positioned at eye, with view vector view, up vector up, and
   // field of view fov.
   void render_aux(int offW, int offH, int side);
   void render(int width, int height, Point3D eye, Vector3D view,
                Vector3D up, float fov, const char* fileName);

   // Add an object into the scene, with material mat.  The function
   // returns a handle to the object node you just added, use the
   // handle to apply transformations to the object.
   SceneDagNode* addObject(SceneObject* obj)
   {
      return addObject(_root, obj);
   }

   // Add an object into the scene with a specific parent node,
   // don't worry about this unless you want to do hierarchical
   // modeling.  You could create nodes with NULL obj and mat,
   // in which case they just represent transformations.
   SceneDagNode* addObject(SceneDagNode* parent, SceneObject* obj);

   // Add a light source.
   LightListNode* addLightSource(LightSource* light);

   // Transformation functions are implemented by right-multiplying
   // the transformation matrix to the node's transformation matrix.

   // Apply rotation about axis 'x', 'y', 'z' angle degrees to node.
   void rotate(SceneDagNode* node, char axis, float angle);

   // Apply translation in the direction of trans to node.
   void translate(SceneDagNode* node, Vector3D trans);

   // Apply scaling about a fixed point origin.
   void scale(SceneDagNode* node, Point3D origin, float factor[3]);
   void traverseEntireScene(Ray3D& ray, bool casting);

   // Flatten scene by pre-applying transforms
   void flattenEntireScene();
   void flattenScene(SceneDagNode* node, Matrix4x4 modelToWorld, Matrix4x4 worldToModel);

   // After intersection, calculate the colour of the ray by shading it
   // with all light sources in the scene.
   void computeShading(Ray3D& ray, int depth, bool getDirectly);

   void printProgress(int percent);

   Colour default_col;

   // Width and height of the viewport.
   int _scrWidth;
   int _scrHeight;
private:
   // Allocates and initializes the pixel buffer for rendering, you
   // could add an interesting background to your scene by modifying
   // this function.
   void initPixelBuffer();

   // Saves the pixel buffer to a file and deletes the buffer.
   void flushPixelBuffer(const char *file_name);

   // Constructs a view to world transformation matrix based on the
   // camera parameters.
   void initInvViewMatrix(Matrix4x4 mat, Point3D eye, Vector3D view,
                           Vector3D up);

   // Traversal code for the scene graph, the ray is transformed into
   // the object space of each node where intersection is performed.
   void traverseScene(SceneDagNode* node, Ray3D& ray, bool casting);

   // Light list and scene graph.
   LightListNode *_lightSource;
   SceneDagNode *_root;

   // Pixel buffer.
   unsigned char* _rbuffer;
   unsigned char* _gbuffer;
   unsigned char* _bbuffer;

public:
   bool soft_shadows;
   bool direct_illumination;
   bool global_illumination;
   bool caustics;
};

#endif // RAYTRACER_H
