#include <math.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <vector>
#include <reducer_opadd.h>

#include "photonmap.h"
#include "config.h"

/* This structure is used only to locate the
 * nearest photons */
struct NearestPhotons
{
	Point3D pos;
	Vector3D normal;
	int max; /* max number of photons we keep track of */
	int found; /* number of photons found so far */
	float max_dist;
	Photon **photons;
};

PhotonMap *createPhotonMap(int max_photons)
{
	PhotonMap *map = new PhotonMap();
	map->stored_photons = 0;
	map->prev_scale = 1;
	map->max_photons = max_photons;

	map->photons = new Photon[max_photons + 1];

	if (map->photons == NULL) {
		fprintf(stderr, "Out of memory initializing photon map\n");
		exit(-1);
	}

	map->bbox_min[0] = map->bbox_min[1] = map->bbox_min[2] = 1e8f;
	map->bbox_max[0] = map->bbox_max[1] = map->bbox_max[2] = -1e8f;

	return map;
}

void destroyPhotonMap(BalancedPhotonMap *map)
{
	delete [] map->photons;
	delete map;
}

/* photonDir returns the direction of a photon
 * at a given surface position */
static void photonDir(Vector3D *dir, const Photon *p)
{
	dir->v[0] = sin(p->theta) * cos(p->phi * 2);
	dir->v[1] = sin(p->theta) * sin(p->phi * 2);
	dir->v[2] = cos(p->theta);
}

/* Store puts a photon into the flat array that will form
 * the final kd-tree.
 *
 * Call this function to store a photon.  */
void storePhoton(
		PhotonMap *map,
		const Colour &power,
		const Point3D &pos,
		const Vector3D &dir)
{
	int i;
	Photon *node;
	if (map->stored_photons >= map->max_photons) {
		Photon *newMap = (Photon *) realloc(map->photons, sizeof(Photon) *
				(2 * map->max_photons + 1));
		if (newMap == NULL) {
			static int done = 0;
			if (!done)
				fprintf(stderr, "Photon Map Full\n");
			done = 1;
			return;
		}
		map->photons = newMap;
		map->max_photons *= 2;
	}

	map->stored_photons++;
	node = &map->photons[map->stored_photons];

	for (i = 0; i < 3; i++) {
		node->pos[i] = pos[i];

		if (node->pos[i] < map->bbox_min[i])
			map->bbox_min[i] = node->pos[i];
		if (node->pos[i] > map->bbox_max[i])
			map->bbox_max[i] = node->pos[i];
	}
	node->power = power;
	node->theta = acos(dir.v[2]);
	node->phi = atan2(dir.v[1], dir.v[0]) / 2;
	node->phi += M_PI * (node->phi < 0);
//	if (node->phi < 0) node->phi += M_PI;
}

/* scalePhotonPower is used to scale the power of all
 * photons once they have been emitted from the light
 * source. scale = 1/(#emitted photons).
 * Call this function after each light source is processed.  */
void scalePhotonPower(PhotonMap *map, const float scale)
{
	cilk_for(int i = map->prev_scale; i <= map->stored_photons; ++i) {
		map->photons[i].power *= scale;
	}
	map->prev_scale = map->stored_photons;
}

// MedianSplit splits the photon array into two separate
// pieces around the median with all photons below the
// the median in the lower half and all photons above
// than the median in the upper half. The comparison
// criteria is the axis (indicated by the axis parameter)
// (inspired by routine in "Algorithms in C++" by Sedgewick)
static void medianSplit(
		Photon **p,
		const int start,               // start of photon block in array
		const int end,                 // end of photon block in array
		const int median,              // desired median number
		const int axis)                // axis to split along
{
#define swap(ph,a,b) { Photon *ph2=ph[a]; ph[a]=ph[b]; ph[b]=ph2; }
	int left = start;
	int right = end;

	while (right > left) {
		const float v = p[right]->pos[axis];
		int i = left - 1;
		int j = right;
		for (;;) {
			while (p[++i]->pos[axis] < v);
			while (p[--j]->pos[axis] > v && j > left);
			if (i >= j)
				break;
			swap(p, i, j);
		}

		swap(p, i, right);
		if (i >= median)
			right = i - 1;
		if (i <= median)
			left = i + 1;
	}
}

// See "Realistic image synthesis using Photon Mapping" chapter 6
// for an explanation of this function
static void balanceSegment(
		PhotonMap *map,
		Photon **pbal,
		Photon **porg,
		const int index,
		const int start,
		const int end)
{
	// Compute new median
	int axis;
	int median = 1;

	while ((4 * median) <= (end - start + 1))
		median += median;

	if ((3 * median) <= (end - start + 1)) {
		median += median;
		median += start - 1;
	} else
		median = end - median + 1;

	// Find axis to split along
	axis = 2;
	if((map->bbox_max[1] + map->bbox_min[2]) >
			(map->bbox_max[2] + map->bbox_min[1])) {
		axis = 1;
	} else if((map->bbox_max[0] + map->bbox_min[1]) >
			(map->bbox_max[1] + map->bbox_min[0]) &&
			(map->bbox_max[0] + map->bbox_min[2]) >
			(map->bbox_max[2] + map->bbox_min[0])) {
		axis = 0;
	}

	// Partition photon block around the median
	medianSplit(porg, start, end, median, axis);

	pbal[index] = porg[median];
	pbal[index]->plane = axis;

	// Recursively balance the left and right block
	if (median > start) {
		// Balance left segment
		if (start +1 < median) {
			const float tmp = map->bbox_max[axis];
			map->bbox_max[axis] = pbal[index]->pos[axis];
			balanceSegment(map, pbal, porg, 2 * index, start, median - 1);
			map->bbox_max[axis] = tmp;
		} else {
			pbal[2 * index] = porg[start];
		}
	}

	if (median < end) {
		// Balance right segment
		if (median + 1 < end) {
			const float tmp = map->bbox_min[axis];
			map->bbox_min[axis] = pbal[index]->pos[axis];
			balanceSegment(map, pbal, porg, 2 * index + 1, median + 1, end);
			map->bbox_min[axis] = tmp;
		} else {
			pbal[2 * index + 1] = porg[end];
		}
	}
}

/* Balance creates a left balanced kd-tree from the flat photon array.
 * This function should be called before the photon map
 * is used for rendering.  */
BalancedPhotonMap *balancePhotonMap(PhotonMap *map)
{
	BalancedPhotonMap *bmap;
	if (map->stored_photons > 1) {
		int i;
		int d, j, foo;
		Photon foo_photon;

		// Allocate two temporary arrays for the balancing procedure
		Photon **pa1 = new Photon*[map->stored_photons + 1];
		Photon **pa2 = new Photon*[map->stored_photons + 1];

		for (i = 0; i <= map->stored_photons; i++)
			pa2[i] = &map->photons[i];

		balanceSegment(map, pa1, pa2, 1, 1, map->stored_photons);
		delete [] pa2;

		// Reorganize balanced kd-tree (make a heap)
		j = 1;
		foo = 1;
		foo_photon = map->photons[j];

		for (i = 1; i <= map->stored_photons; i++) {
			d = pa1[j] - map->photons;
			pa1[j] = NULL;
			if (d != foo)
				map->photons[j] = map->photons[d];
			else {
				map->photons[j] = foo_photon;

				if (i < map->stored_photons) {
					for (;foo <= map->stored_photons; foo++)
						if (pa1[foo] != NULL)
							break;
					foo_photon = map->photons[foo];
					j = foo;
				}
				continue;
			}
			j = d;
		}
		delete [] pa1;
	}

	bmap = new BalancedPhotonMap();
	bmap->stored_photons = map->stored_photons;
	bmap->half_stored_photons = map->stored_photons / 2 - 1;
	bmap->photons = map->photons;
	delete map;

	return bmap;
}

/* Locate_photons finds the nearest photons in the
 * photon map given the parameters in np */
static bool locatePhotons(
		BalancedPhotonMap *map,
		NearestPhotons *const np,
		const int index,
		const Vector3D &normal)
{
	Photon *p = &map->photons[index];
	if (index < map->half_stored_photons)
	{
		float dist1 = np->pos[p->plane] - p->pos[p->plane];

		if (dist1 > 0.0) {
			// If dist1 is positive search right plane
			if(!locatePhotons(map, np, 2 * index + 1, normal))
				return false;
			if(dist1 * dist1 < np->max_dist && !locatePhotons(map, np, 2 * index , normal))
					return false;
		} else {
			// Else, dist1 is negative search left first
			if(!locatePhotons(map, np, 2 * index , normal))
				return false;
			if(dist1 * dist1 < np->max_dist && !locatePhotons(map, np, 2 * index + 1 , normal))
				return false;
		}
	}

	// Compute squared distance between current photon and np->pos
	Vector3D del = p->pos - np->pos;
	float dist2 = del.dot(del);

	// Adjust the distance for photons that are not on the same plane as this
	// point.
	float discFix = normal.dot(del);
	discFix = fabs(discFix);
	dist2 += discFix * dist2 * 0.010;

	if (dist2 < np->max_dist) {
		// We found a photon :) Check length and insert
		if(np->found < np->max) {
			np->photons[np->found++] = p;
		} else {
			return false;
		}
	}

	return true;
}

/* Irradiance_estimate computes an irradiance estimate
 * at a given surface position */
void irradianceEstimate(
		BalancedPhotonMap *map,
		Colour *irrad,                // Returned irradiance
		const Point3D &pos,           // Surface position
		const Vector3D &normal,       // Surface normal at pos
		const float max_dist,        // Max distance to look for photons
		const int nphotons)           // Number of photons to use
{
	NearestPhotons np;

	// We used to use allocas here, but cilk doesn't like those, so we use
	// stack-allocated vectors which automatically free their storage on return.
	np.pos = pos;
	np.max = nphotons;
	np.max_dist = max_dist * max_dist;
	np.found = 0;
	np.photons = new Photon*[nphotons];//new max_photon_heap(nphotons);

	// Locate the nearest photons
	locatePhotons(map, &np, 1 , normal);

	// If less than 2 photons return
	/* why 2? */
	if (np.found < 20) {
		return;
	}

	// Sum irradiance from all photons
	/* in order to cilk_for this, we need a hyperobject wrapper around the color
	 * otherwise we have race conditions */
	cilk::reducer_opadd<float> r;
	cilk::reducer_opadd<float> g;
	cilk::reducer_opadd<float> b;
	cilk_for(int i = 0; i < np.found; i++) {
		Vector3D pdir;
		const Photon *p = np.photons[i];
		photonDir(&pdir, p);
		if (pdir.dot(normal) < 0.0) {
			r += p->power.r;
			g += p->power.g;
			b += p->power.b;
		}
	}
	*irrad = Colour(r.get_value(), g.get_value(), b.get_value());

	// Take into account (estimate of) density
	/* 1.0f / M_PI = 0.3183098861837906715 */
	*irrad *= 0.3183098861837906715 / (np.max_dist);
	delete np.photons;
}
