#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#include "util.h"
#include "config.h"

/* This is the photon.  This should be as small as possible to minimize the size
 * of the kd-tree (for spatial locality). */
struct Photon
{
	Point3D  pos;                // Photon position
	Colour power;                // Photon power
	int plane;                   // Splitting plane for kd-tree (one of three
								 // values: 0, 1, or 2)
	float theta;
	float phi;                  // Since we need only very little accuracy for
								// the photon direction, we can save some
								// space by storing it in spherical coords.
	float dist2;				// Used for nearestphoton calculations
};

struct BalancedPhotonMap
{
   Photon *photons;       // Array of photons (making up the kd-tree)
   int half_stored_photons;
   int stored_photons;
};


/* The actual photon map structure */
struct PhotonMap
{
   float bbox_min[3];          // use bbox_min;
   float bbox_max[3];          // use bbox_max;

   int stored_photons;
   Photon *photons;       // Array of photons (making up the kd-tree)
   int half_stored_photons;

   // The photon map MUST be the same as the balanced map up to this point.
   // What follows is only used when the map is created...
   int max_photons;
   int prev_scale;
};


PhotonMap *createPhotonMap(int max_photons);

void storePhoton(PhotonMap *map,
                 const Colour &power,          // Photon power
                 const Point3D &pos,           // Photon position
                 const Vector3D &dir);         // Photon direction

void scalePhotonPower(PhotonMap *map,
                      const float scale);   // 1 / (number of emitted photons)

BalancedPhotonMap *balancePhotonMap(PhotonMap *map);  // Balance the kd-tree

void irradianceEstimate(BalancedPhotonMap *map,
                        Colour *irrad,          // Returned irradiance
                        const Point3D &pos,     // Surface position
                        const Vector3D &normal, // Surface normal at pos
                        const float max_dist,  // Max dist to look for photons
                        const int nphotons );   // Number of photons to use

void destroyPhotonMap(BalancedPhotonMap *map);

#endif // PHOTONMAP_H
