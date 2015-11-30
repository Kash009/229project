#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

// Example representation of an autofocus zone.
class AfZone {
	public:
	  // from the data file
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};
class Lens {
   public:
      bool intersect(Ray &ray, float nextmaterial);
      float radius;
      float z;
      float refraction;
      float haperture;

};
class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
	  const string &autofocusfile,
      float filmdiag,
	  Film *film);
   ~RealisticCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
   void UpdateDistance(float dist);
   void  ParseLens(const string& filename, float apertured);
   void  ParseAfZones(const string& filename);

private:
   bool  autofocus;
   vector<AfZone> afZones;
   vector<Lens> lensSet;
   float filmSizeConst;
   float diskradius;
   float weightConst;
   float filmDist;
   float F[33];
   //float threshold = 0.9;
   //int nsamples = 9;
   float filmDistStart;
   float filmDistEnd;
   int xres;
   int yres;
   float ShutterOpen;
   float ShutterClose;
   Film * film;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
