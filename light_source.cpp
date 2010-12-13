#include <cmath>
#include <iostream>
#include <assert.h>
#include "light_source.h"
#include "photonmap.h"
#include "raytracer.h"
#include "random.h"
#include "config.h"

SquarePhotonLight::SquarePhotonLight(Colour col, Raytracer *raytracer ) :
	col(col), raytracer(raytracer), icache(ICACHE_TOLERANCE,
			ICACHE_MIN_SPACING)
{
	light_col = col;
}

SquarePhotonLight::~SquarePhotonLight() {
	destroyPhotonMap(bmap);
	destroyPhotonMap(cmap);
}

/* seems like this is where Matrix4x4 is used */
void SquarePhotonLight::initTransformMatrix(Matrix4x4 mat, Vector3D& w)
{
	Vector3D u, v;

	if ((fabs(w.v[0]) < fabs(w.v[1])) && (fabs(w.v[0]) < fabs(w.v[2]))) {
		v.v[0] = 0;
		v.v[1] = w.v[2];
		v.v[2] = -w.v[1];
	} else if (fabs(w.v[1]) < fabs(w.v[2])) {
		v.v[0] = w.v[2];
		v.v[1] = 0;
		v.v[2] = -w.v[0];
	} else {
		v.v[0] = w.v[1];
		v.v[1] = -w.v[0];
		v.v[2] = 0;
	}
	v.normalize();
	u = v.cross(w);

	mat[0][0] = u.v[0];
	mat[1][0] = u.v[1];
	mat[2][0] = u.v[2];
	mat[0][1] = v.v[0];
	mat[1][1] = v.v[1];
	mat[2][1] = v.v[2];
	mat[0][2] = w.v[0];
	mat[1][2] = w.v[1];
	mat[2][2] = w.v[2];
}

void SquarePhotonLight::globalIllumination(Ray3D& ray, bool getDirectly)
{
	// If we're not already in the irradiance cache, compute the irradiance via
	// monte carlo methods.
	if (!icache.getIrradiance(ray.intersection.point, ray.intersection.normal, &ray.col)) {
		Colour c;

		ray.col = Colour(0, 0, 0);

		int N = MONTE_CARLO_STRATIFICATION_N;
		int M = MONTE_CARLO_STRATIFICATION_M;
		int hits = 0;
		float r0 = 0;

		Matrix4x4 basis;
		initTransformMatrix(basis, ray.intersection.normal);

		// Stratification
		for (int i = 0; i < N; i++) {
			float phi = 2 * M_PI * ((float) i + r.Random1()) / N;
			float sinPhi = sin(phi);
			float cosPhi = cos(phi);

			for (int j = 0; j < M; j++) {
				float cosTheta = rsqrtss_sqrt(1 - (((float) j + r.Random1()) / M));
				float theta = acos(cosTheta);
				float sinTheta = sin(theta);

				Vector3D v = Vector3D(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
				v = v.transform(basis);
				v.normalize();

				Ray3D new_ray = Ray3D(ray.intersection.point, v);
				raytracer->traverseEntireScene(new_ray, true);

				if (!new_ray.intersection.none) {
					raytracer->computeShading(new_ray, 3, true);
					ray.col += new_ray.col;

					r0 += 1 / new_ray.intersection.t_value;
					hits++;
				}
			}
		}

		ray.col *= 1.0 / hits;
		r0 = 1 / r0;

		if (hits == N * M) {
			icache.insert(ray.intersection.point, ray.intersection.normal, r0, ray.col);
		}
	}
	ray.col *= ray.intersection.mat->diffuse;
}

void SquarePhotonLight::causticIllumination(Ray3D& ray)
{
	// Caustics
	Colour caus_col;

	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	if (raytracer->soft_shadows) {
		irradianceEstimate(cmap,
				&caus_col,
				ray.intersection.point,
				normal,
				CAUSTICS_SOFT_MAX_DISTANCE,
				CAUSTIC_SOFT_MAX_PHOTONS);
	} else {
		irradianceEstimate(cmap,
				&caus_col,
				ray.intersection.point,
				normal,
				CAUSTICS_MAX_DISTANCE,
				CAUSTIC_MAX_PHOTONS);
	}
	float cosTheta12 = rsqrtss_sqrt(-(ray.dir.dot(ray.intersection.normal)));
	caus_col *= ray.intersection.mat->diffuse;

	ray.col += caus_col * cosTheta12;
}

void SquarePhotonLight::directIllumination(Ray3D& ray)
{
	// Direct illumination
	int N = 1;
	int M = 1;

	if (raytracer->soft_shadows) {
		N = NUM_SOFT_SHADOW_RAYS_IN_EACH_DIM;
		M = NUM_SOFT_SHADOW_RAYS_IN_EACH_DIM;
	}

	Colour direct_col = Colour(0, 0, 0);
	float dx = 1.0 / (N + 1);
	float dz = 1.0 / (M + 1);

	// Loop for soft shadows
	for (int i = 1; i <= N; i++) {
		float x;
		if (raytracer->soft_shadows) {
			float rand = r.Random1();
			float dx_rand = rand - floor(rand);
			x = ((i + dx_rand) * dx) * 30 - 15;
		} else {
			x = i * dx * 30 - 15;
		}

		for (int j = 1; j <= M; j++) {
			float z;
			if (raytracer->soft_shadows) {
				float rand = r.Random1();
				float dz_rand = rand - floor(rand);
				z = (((float) j + dz_rand) * dz) * 30 - 15;
			} else {
				z = i * dz * 30 - 15;
			}

			Vector3D L = Point3D(0 + x, 50, 0 + z) - ray.intersection.point;
			float l = rsqrtss_sqrt(L.v[0] * L.v[0] + L.v[1] * L.v[1] + L.v[2] * L.v[2]);
			L.normalize();

			Vector3D LN = Vector3D(0, -1, 0);

			float scale = -LN.dot(L);
			scale *= (scale >= 0);
//			scale = scale < 0 ? 0 : scale;
			scale /= l * l * 1.5 * M_PI;

			Ray3D new_ray = Ray3D(ray.intersection.point, L);
			raytracer->traverseEntireScene(new_ray, false);

			if (!new_ray.intersection.none && new_ray.intersection.mat->light) {
				Vector3D R = 2.0 * ray.intersection.normal.dot(L) * ray.intersection.normal - L;

				float NdotL = L.dot(ray.intersection.normal);
				float RdotV = -(R.dot(ray.dir));
				NdotL *= (NdotL >= 0);
//				NdotL = NdotL < 0 ? 0 : NdotL;
				RdotV *= (RdotV >= 0);
//				RdotV = RdotV < 0 ? 0 : RdotV;
				RdotV *= ray.dir.dot(ray.intersection.normal) <= 0;
				/*
				if (ray.dir.dot(ray.intersection.normal) > 0) {
					RdotV = 0;
				}
				*/
				direct_col += (light_col * scale * (ray.intersection.mat->diffuse * NdotL + ray.intersection.mat->specular * pow(RdotV, ray.intersection.mat->specular_exp)));
			}
		}
	}
	direct_col = direct_col  / (N * M);
	ray.col += direct_col;
}

void SquarePhotonLight::shade(Ray3D& ray, bool getDirectly)
{
	// Don't bother shading lights, just give it the color of the light source.
	if (ray.intersection.mat->light) {
		ray.col = ray.intersection.mat->diffuse;
		return;
	}

	// getDirectly means we visualize the photon map directly.
	if (getDirectly) {

		Vector3D normal = ray.intersection.normal;
		normal.normalize();

		irradianceEstimate(bmap,
				&ray.col,
				ray.intersection.point,
				normal,
				INDIRECT_MAX_DISTANCE,
				INDIRECT_MAX_PHOTONS);
		ray.col = ray.col * ray.intersection.mat->diffuse;
		return;
	}

	// Only look at objects in the photon map
	if (ray.intersection.mat->isDiffuse) {

		if (raytracer->global_illumination)
			globalIllumination(ray, getDirectly);

		if (raytracer->caustics)
			causticIllumination(ray);
	}

	if (raytracer->direct_illumination)
		directIllumination(ray);
}

Vector3D SquarePhotonLight::getRandLambertianDir(Vector3D& normal)
{
	Vector3D v;

	float phi = 2 * M_PI * r.Random1();
	float sinPhi = sin(phi);
	float cosPhi = cos(phi);

	float cosTheta = rsqrtss_sqrt(r.Random1());
	float theta = acos(cosTheta);
	float sinTheta = sin(theta);

	Matrix4x4 basis;
	initTransformMatrix(basis, normal);

	v = Vector3D(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
	v = v.transform(basis);
	v.normalize();

	return v;
}

void SquarePhotonLight::tracePhotons(int num, int caustics_num)
{
	// Flatten scene before render
	raytracer->flattenEntireScene();

	PhotonMap *map = createPhotonMap(num * 4);

	cout << "Emitting global illumination photons: 0%% ";
	int i = 0;
	while (i < num) {
		// Calculate the start pos and direction of the photon
		Point3D p = Point3D(r.Random2() * 16, 49.99, r.Random2() * 16);
		Vector3D v = Vector3D(0, -1, 0);
		v = getRandLambertianDir(v);
		Ray3D ray = Ray3D(p, v);
		ray.col = col;

		int count = 0;
		while (ray.col.max() > 0.1 && count++ < 100) {
			raytracer->traverseEntireScene(ray, true);

			if (ray.intersection.none) {
				i--;
				break;
			}

			Vector3D dir_norm = ray.dir;
			dir_norm.normalize();
			if (ray.intersection.mat->isDiffuse) {
				storePhoton(map,
						ray.col,
						ray.intersection.point,
						dir_norm);
			}

			float ran = r.Random1();

			Colour c = ray.col * ray.intersection.mat->diffuse;
			float P = c.max() / ray.col.max();
			if (ran < P) {
				// Diffuse reflection
				v = getRandLambertianDir(ray.intersection.normal);
			} else {
				ran -= P;
				c = ray.col * ray.intersection.mat->specular;
				P = c.max() / ray.col.max();
				if (ran < P) {
					// Specular reflection
					v = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
						ray.intersection.normal;
				} else {
					ran -= P;
					c = ray.col * ray.intersection.mat->refractive;
					P = c.max() / ray.col.max();
					if (ran < P) {
						// Refraction
						float n;

						if (ray.dir.dot(ray.intersection.normal) < 0) {
							n = 1 / ray.intersection.mat->refr_index;
						} else {
							ray.intersection.normal = -ray.intersection.normal;
							n = ray.intersection.mat->refr_index;
						}

						float cosI = ray.intersection.normal.dot(ray.dir);
						float sinT2 = n * n * (1.0 - cosI * cosI);

						if (sinT2 < 1.0) {
							v = n * ray.dir - (n * cosI + rsqrtss_sqrt(1.0 - sinT2)) *
								ray.intersection.normal;
						} else {
							// Total internal reflection
							v = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
								ray.intersection.normal;
						}
					} else {
						// Absorption
						break;
					}

				}
			}
			ray.origin = ray.intersection.point;
			ray.dir = v;
			ray.col = c / P;
			ray.intersection.none = true;
			ray.intersection.t_value = FLT_MAX;
		}
		raytracer->printProgress((int) (((i + 1) * 100.0) / num));
		i++;
	}
	scalePhotonPower(map, 1.0 / i);
	bmap = balancePhotonMap(map);

	// Caustics
	map = createPhotonMap(caustics_num);

	int emitted = 0;
	i = 0;
	cout << endl << "Emitting caustics illumination photons: 0%% ";

	while (i < caustics_num) {
		// Calculate the start pos and direction of the photon
		Point3D p;

		if (raytracer->soft_shadows) {
			p = Point3D(r.Random2() * 18, 49.99, r.Random2() * 18);
		} else {
			p = Point3D(0, 49.99, 0);
		}
		Vector3D v = Vector3D(0, -1, 0);
		v = getRandLambertianDir(v);
		Ray3D ray = Ray3D(p, v);
		ray.col = col;
		emitted++;

		raytracer->traverseEntireScene(ray, true);

		if (!ray.intersection.none && ray.intersection.mat->isSpecular) {
			while (1) {
				raytracer->traverseEntireScene(ray, true);

				if (ray.intersection.none) break;

				Vector3D dir_norm = ray.dir;
				dir_norm.normalize();
				if (ray.intersection.mat->isDiffuse) {
					storePhoton(map,
							ray.col,
							ray.intersection.point,
							dir_norm);

					i++;
					break;
				}

				float ran = r.Random1();
				Colour c = ray.col * ray.intersection.mat->specular;
				float P = c.max() / ray.col.max();

				// Specular reflection
				if (ran < P) {

					ray.intersection.normal = -ray.intersection.normal;

					v = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
						ray.intersection.normal;

				} else {
					ran -= P;
					c = ray.col * ray.intersection.mat->refractive;
					P = c.max() / ray.col.max();
					if (ran < P) {
						// Refraction
						float n;

						if (ray.dir.dot(ray.intersection.normal) < 0) {
							n = 1 / ray.intersection.mat->refr_index;
						} else {
							ray.intersection.normal = -ray.intersection.normal;
							n = ray.intersection.mat->refr_index;
						}

						float cosI = ray.intersection.normal.dot(ray.dir);
						float sinT2 = n * n * (1.0 - cosI * cosI);

						if (sinT2 < 1.0) {
							v = n * ray.dir - (n * cosI + rsqrtss_sqrt(1.0 - sinT2)) *
								ray.intersection.normal;
						} else { // Total internal reflection.
							v = ray.dir - (2 * ray.dir.dot(ray.intersection.normal)) *
								ray.intersection.normal;
						}
					} else {
						// Absorption
						i++;
						break;
					}
				}
				ray.origin = ray.intersection.point;
				ray.dir = v;
				ray.col = c / P;
				ray.intersection.none = true;
				ray.intersection.t_value = FLT_MAX;
			}

			raytracer->printProgress((int) (((i + 1) * 100.0) / caustics_num));
		}
	}
	cout << endl;
	scalePhotonPower(map, 1.0 / emitted);
	cmap = balancePhotonMap(map);
}
