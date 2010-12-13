#include "raytracer.h"
#include <stdio.h>

#ifdef __cilkplusplus
#define MAIN cilk_main
#else
#define MAIN main
#endif

// Constructs the scene to be rendered, casts photons, makes the irradiance
// map, and renders the scene.
int MAIN(int argc, char* argv[])
{
   // Rendering options
   int scene = 1;
   int size = 600;
   const char *outfile = "output.bmp";
   bool direct_illumination = 1;
   bool global_illumination = 1;
   bool caustics = 1;

   // Parse arguments
   int c;
   while ((c = getopt(argc, argv, "s:n:o:d:g:c:h")) != -1)
      switch (c) {
      case 's':
         scene = atoi(optarg);
         break;
      case 'n':
         size = atoi(optarg);
         break;
      case 'o':
         outfile = optarg;
         break;
      case 'd':
         if (!strcmp(optarg, "off") || !strcmp(optarg, "0"))
            direct_illumination = 0;
         else
            direct_illumination = 1;
         break;
      case 'g':
         if (!strcmp(optarg, "off") || !strcmp(optarg, "0"))
            global_illumination = 0;
         else
            global_illumination = 1;
         break;
      case 'c':
         if (!strcmp(optarg, "off") || !strcmp(optarg, "0"))
            caustics = 0;
         else
            caustics = 1;
         break;
      case 'h':
      case '?':
         printf("Usage: raytracer [options]\n"
                "  -s <scene>     Scene to render, default 1\n"
                "  -n <size>      Image pixel width/height, default 600\n"
                "  -o <outfile>   Output filename\n"
                "  -d (on|off)    Enable/disable direct illumination\n"
                "  -g (on|off)    Enable/disable global illumination\n"
                "  -c (on|off)    Enable/disable caustic illumination\n"
                "  -h             Print help/usage\n");
         return 1;
      default:
         abort();
      }

   printf("Scene: %d\n", scene);
   printf("Size: %d\n", size);
   printf("Output file: %s\n", outfile);
   printf("Direct illumination: %s\n", (direct_illumination ? "on" : "off"));
   printf("Global illumination: %s\n", (global_illumination ? "on" : "off"));
   printf("Caustics: %s\n", (caustics ? "on" : "off"));

   // Camera parameters.
   Point3D eye(0, 20, 100);
   Vector3D gaze(0, 0, -1);
   Vector3D up(0, 1, 0);
   float fov = 20;

   // Define materials for shading.
   Material red(Colour(0.92941176,  0.52156863,  0.14117647),
                Colour(0.0, 0.0, 0.0), 4, Colour(0.0, 0.0, 0.0), 1, true);
   Material green(Colour(0, 0.8, 0), Colour(0.0, 0.0, 0.0), 4,
                  Colour(0.0, 0.0, 0.0), 1, true);
   Material blue(Colour(0.14117647,  0.52156863,  0.92941176),
                 Colour(0.0, 0.0, 0.0), 4, Colour(0.0, 0.0, 0.0), 1, true);
   Material grey(Colour(0.7, 0.7, 0.7), Colour(0.0, 0.0, 0.0), 4,
                 Colour(0.0, 0.0, 0.0), 1, true);
   Material yellow(Colour(0xEE/255.,  0xF7/255.,  0x45/255.),
                   Colour(0.0, 0.0, 0.0), 4, Colour(0.0, 0.0, 0.0), 1, true);
   Material magenta(Colour(0xE0/255.,  0x1B/255.,  0x4C/255.),
                    Colour(0.0, 0.0, 0.0), 4, Colour(0.0, 0.0, 0.0), 1, true);

   Material glass(Colour(0.0, 0.0, 0.0), Colour(0.07, 0.07, 0.07), 5,
                  Colour(0.89, 0.89, 0.89), 1.8, false, true);
   Material water(Colour(0.0, 0.00, 0.0), Colour(0.07, 0.07, 0.07), 5,
                  Colour(0.93, 0.93, 0.93), 1.3, false, true);
   Material mirror(Colour(0.0, 0.0, 0.0), Colour(0.7, 0.7, 0.7), 25,
                   Colour(0.0, 0.0, 0.0), 1.7, false, true);
   Material lightmat(Colour(1.0, 1.0, 1.0), Colour(0.0, 0.0, 0.0), 55,
                     Colour(0.0, 0.0, 0.0), 1.0, false, true, true);

   if (scene == 1) {
      Raytracer raytracer(true, direct_illumination, global_illumination, caustics);

      // Add objects to the scene
      raytracer.addObject(new Cube(IN, &grey, &grey, &grey, &grey, &blue, &red));
      raytracer.addObject(new Square(&lightmat));

      // Mirror surface sphere
      SceneDagNode* sphere1 = raytracer.addObject(new Sphere(&mirror));
      raytracer.translate(sphere1, Vector3D(-27, -30, -30));
      // Solid glass sphere
      SceneDagNode* sphere2 = raytracer.addObject(new Sphere(&glass));
      raytracer.translate(sphere2, Vector3D(29, -30, -5));

      // The light source at the top
      SquarePhotonLight *light =
          new SquarePhotonLight(Colour(15000.0, 15000.0, 15000.0), &raytracer);

      light->tracePhotons(NUM_PHOTONS, SCENE_1_NUM_CAUSTIC_PHOTONS);
      raytracer.addLightSource(light);

      raytracer.render(size, size, eye, gaze, up, fov, outfile);
   } else if (scene == 2) {
      Raytracer raytracer(false, direct_illumination, global_illumination, caustics);

      // Add objects to the scene
      raytracer.addObject(new Cube(IN, &grey, &grey, &grey, &grey, &green, &blue));
      raytracer.addObject(new Square(&lightmat));

      // Mirror surface sphere
      SceneDagNode* sphere1 = raytracer.addObject(new Sphere(&mirror));
      raytracer.translate(sphere1, Vector3D(-27, -30, -30));

      // Water surface
      raytracer.addObject(new DisplacedSurface(&water, 22.0));

      // The light source at the top
      SquarePhotonLight *light =
          new SquarePhotonLight(Colour(15000.0, 15000.0, 15000.0), &raytracer);

      light->tracePhotons(NUM_PHOTONS, SCENE_2_NUM_CAUSTIC_PHOTONS);
      raytracer.addLightSource(light);

      raytracer.render(size, size, eye, gaze, up, fov, outfile);
   } else if (scene == 3) {
      Raytracer raytracer(false, direct_illumination, global_illumination, caustics);

      // Add objects to the scene
      raytracer.addObject(new Cube(IN, &grey, &grey, &grey, &grey, &yellow, &magenta));
      raytracer.addObject(new Square(&lightmat));

      // Mirror surface sphere
      SceneDagNode* sphere1 = raytracer.addObject(new Sphere(&mirror));
      raytracer.translate(sphere1, Vector3D(27, 30, 30));
      // Solid glass sphere
      SceneDagNode* sphere2 = raytracer.addObject(new Sphere(&glass));
      raytracer.translate(sphere2, Vector3D(-29, 30, 5));

      // Water surface
      raytracer.addObject(new DisplacedSurface(&water, -22.0));

      // The light source at the top
      SquarePhotonLight *light =
            new SquarePhotonLight(Colour(15000.0, 15000.0, 15000.0), &raytracer);

      light->tracePhotons(NUM_PHOTONS, SCENE_2_NUM_CAUSTIC_PHOTONS);
      raytracer.addLightSource(light);

      raytracer.render(size, size, eye, gaze, up, fov, outfile);
   } else {
      std::cout << "Bad scene number.\n";
      return 1;
   }

   return 0;
}
