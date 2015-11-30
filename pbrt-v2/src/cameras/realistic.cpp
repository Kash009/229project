// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "cameras/vdb.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>

using namespace std;
RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {	   
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, autofocusfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   filmDist(filmdistance),
								   film(f)
{
	ParseLens(specfile, aperture_diameter_);
	UpdateDistance(filmDist);
	
	autofocus = false;
	if (autofocusfile.compare("") != 0){
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
	xres = f->xResolution;
	yres = f->yResolution;
	filmSizeConst = filmdiag/sqrtf(xres*xres + yres*yres);
	filmDistStart = 0.2*filmDist;
	filmDistEnd = 2.4*filmDist;
	//printf("filmDist(%f),start(%f),end(%f)\n", filmDist, filmDistStart, filmDistEnd);
}

// parses the AF zone file
void RealisticCamera::ParseAfZones(const string& filename)
{
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		afZones.resize(afZones.size()+1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }

	printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

// parses the lens specs file

void RealisticCamera::ParseLens(const string& filename, float apertured)
{
  ifstream lensfile(filename.c_str());
   if (!lensfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }
   char line[512];
   while (!lensfile.eof()) {
      lensfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		lensSet.resize(lensSet.size()+1);
		Lens& lens = lensSet[lensSet.size()-1];
		sscanf(line, "%f %f %f %f\n", &lens.radius, &lens.z, &lens.refraction, &lens.haperture);
		lens.haperture*=0.5;

		if(lens.radius < 0.0001 && -0.0001 < lens.radius){
			lens.haperture = 0.5*apertured; 
		}
	   }
	}

	//set distance to film as first lens' offset
	reverse(lensSet.begin(),lensSet.end());
	Lens &lens = lensSet[0];
	diskradius = lens.haperture;
	lens.z=filmDist;
	printf("Read in %zu lens specs from %s\n", lensSet.size(), filename.c_str());
}


//Before this function, ray is positioned at previous lens, origin is center of previous lens.
bool inAperture(float x, float y, float r){
		//flower
		//float theta = acosf(x/(x*x+y*y));
		//return x*x +y*y < 0.36*(powf(r*(0.4+0.6*pow(1-fabsf(cos(2.5f*theta)),0.75)),2));
		return x*x+y*y < r*r;

}


bool Lens::intersect(Ray& ray, float inext){
	
	//get intersection with film plane side of lens
	Point intersection;
	
	//if this is stop, directly check intersection without spherical constraints
	if(radius < 0.0001 && -0.0001 < radius){
		float time = (z-(ray.o).z)/((ray.d).z);
		intersection = ray(time);
		if(!inAperture(intersection.x,intersection.y, haperture))return false;	
	}else{
		//find intersection with sphere
		Point center(0.f,0.f,-radius+z);
		Vector tocenter(center.x-(ray.o).x,center.y-(ray.o).y, center.z -(ray.o).z);
		float dintersect = Dot(ray.d,tocenter);
		float check = dintersect*dintersect + radius*radius - tocenter.LengthSquared();
		if(check <0)return false;

		//find smallest positive time
		float time = dintersect - sqrtf(check);
		if(time <0)time += 2*sqrtf(check);
		if(time <0)return false;
		intersection = ray(time);

		//check if is within aperture
		if((intersection.x)*(intersection.x)+(intersection.y)*(intersection.y) > haperture*haperture) return false;

		//calculate normal of point on lens
		Vector n(intersection);
		if(radius >0)n *= -1;
		n.z = -1*sqrtf(radius*radius-n.x*n.x-n.y*n.y);
		n = Normalize(n);

		//Refract
		if(inext < 0.0001)inext = 1.0f;
		float mu = refraction/inext;
		float dotp = Dot(n,ray.d);
		float coesq = 1-mu*mu*(1-dotp*dotp);
		if(coesq <0) return false;	//internal reflection
		float ncoeff = -mu*dotp - sqrtf(coesq);
		ray.d = mu* ray.d + ncoeff * n;
	}
	//move ray origin to proper point aka origin is at z-intercept of this lens
	ray.o = intersection;
	ray.o.z -= z;
	ray = Ray(ray.o, ray.d, 0.0f);
	return true;
}

RealisticCamera::~RealisticCamera(){}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{		
     	//sample film point
       float x = filmSizeConst*(sample.imageX - 0.5*xres);
       float y = filmSizeConst*(sample.imageY - 0.5*yres);
       Point filmp(-x,y,0);

		//sample disk      
       float u,v;
       ConcentricSampleDisk(sample.lensU, sample.lensV, &u, &v);
       u *= diskradius;
       v *= diskradius;
       Point diskp(u,v,filmDist);
      
       Vector dir(diskp.x - filmp.x, diskp.y -filmp.y, diskp.z - filmp.z);
       dir = Normalize(dir);
     //  float weight = 0.5f*(1-((x*x+y*y+filmDist*filmDist-diskradius*diskradius)/sqrtf(powf(x*x+y*y+filmDist*filmDist+u*u+v*v,2)-4*diskradius*diskradius*(x*x+y*y))));

       Ray r(filmp, dir, 0.0f);


       //start updating ray based on refractions
       int size = lensSet.size();
       for(int i = 0; i < size-1; i++){
       		Lens l = lensSet[i];
       		float inext = (lensSet[i+1]).refraction;
       		if(!(l.intersect(r,inext)))
       			return 0.f;
       }

       //last ray update into world, air coeff 1
       Lens final = lensSet[size-1];
       if(!(final.intersect(r,1.0f))) return 0.f;
    CameraToWorld(r,ray);
	ray->d = Normalize(ray->d);

	//return weight;
  return weightConst*pow(dir.z,4);
}

void RealisticCamera::UpdateDistance(float dist){
	Lens &first = lensSet[0];
	if(first.radius >0){
		float delta = fabsf(first.radius) - sqrtf(first.radius*first.radius - first.haperture*first.haperture);
		diskradius *= filmDist/(filmDist-delta);
	}
	weightConst = (diskradius*diskradius*M_PI)/((filmDist*filmDist));
	filmDist = dist;
	first.z=filmDist;
}

void  RealisticCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {
	// 1. Modify this code so that it can adjust film plane of the camera
	// 2. Use the results of raytracing to evaluate whether the image is in focus
	// 3. Search over the space of film planes to find the best-focused plane.
	
	//UpdateDistance(62.0f);
	//return;
	if(!autofocus)
		return;
	int slice = 32;
	vector<float> filmdists;
	filmdists.assign(afZones.size(),3*filmDist);
	for(int recurse = 0; recurse <2; recurse++){
		float increment = (1.0f/(slice+1))*(filmDistEnd - filmDistStart);

		for (size_t i=0; i<afZones.size(); i++) {
			AfZone & zone = afZones[i];
			for(int j = 0; j < slice+1; j++){
				printf("start%d-",j);
				UpdateDistance(filmDistStart+j*increment);
				RNG rng;
				MemoryArena arena;
				Filter * filter = new BoxFilter(.5f,.5f);
				const float crop[] = {zone.left,zone.right,zone.top,zone.bottom};
				ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,"foo.exr",false);
				int xstart,xend,ystart,yend;
				sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);

				StratifiedSampler sampler(xstart, xend, ystart, yend,
				                          16, 16, true, ShutterOpen, ShutterClose);
				// Allocate space for samples and intersections
				int maxSamples = sampler.MaximumSampleCount();
				Sample *samples = origSample->Duplicate(maxSamples);
				RayDifferential *rays = new RayDifferential[maxSamples];
				Spectrum *Ls = new Spectrum[maxSamples];
				Spectrum *Ts = new Spectrum[maxSamples];
				Intersection *isects = new Intersection[maxSamples];
				
					// Get samples from _Sampler_ and update image
				int sampleCount;
				weightConst = 1;
				
				while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
					// Generate camera rays and compute radiance along rays
					for (int i = 0; i < sampleCount; ++i) {
						// Find camera ray for _sample[i]_

						float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
						rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));


						// Evaluate radiance along camera ray

						if (rayWeight > 0.f)
							Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
															 arena, &isects[i], &Ts[i]);
						else {
							Ls[i] = 0.f;
							Ts[i] = 1.f;
						}

						// Issue warning if unexpected radiance value returned
						if (Ls[i].HasNaNs()) {
							Error("Not-a-number radiance value returned "
								  "for image sample.  Setting to black.");
							Ls[i] = Spectrum(0.f);
						}
						else if (Ls[i].y() < -1e-5) {
							Error("Negative luminance value, %f, returned"
								  "for image sample.  Setting to black.", Ls[i].y());
							Ls[i] = Spectrum(0.f);
						}
						else if (isinf(Ls[i].y())) {
							Error("Infinite luminance value returned"
								  "for image sample.  Setting to black.");
							Ls[i] = Spectrum(0.f);
						}

					}

					// Report sample results to _Sampler_, add contributions to image
					if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
					{
						for (int i = 0; i < sampleCount; ++i)
						{

							sensor.AddSample(samples[i], Ls[i]);

						}
					}
					// Free _MemoryArena_ memory from computing image sample values
					arena.FreeAll();
				}

				float * rgb;
				int width;
				int height;
				sensor.WriteRGB(&rgb,&width,&height,1.f);
				float f=0;
				float averageI=0;

				for(int y = 1; y <height-1; y++){
					for(int x = 1; x <width-1; x++){
						int offset = 3*(y*width + x);
						int yoff = 3*width;
						averageI += 0.2989*rgb[offset]+0.5870*rgb[offset+1]+0.1140*rgb[offset+2];
						/*directly using sobel operator (convolution with
						| .5  1  .5|
						|  0  0   0|
						|-.5 -1 -.5| and its transpose. catches horizontal/vertical
						)*/

						float mlx = fabsf(0.2989*(rgb[offset+3]
								-rgb[offset-3]
								
								+0.5*(rgb[offset+yoff+3]
								-rgb[offset+yoff-3])
								
								+0.5*(rgb[offset-yoff+3]
								-rgb[offset-yoff-3])));

						mlx += fabsf(0.5870*(rgb[offset+3+1]
								-rgb[offset-3+1]
								
								+0.5*(rgb[offset+yoff+3+1]
								-rgb[offset+yoff-3+1])
								
								+0.5*(rgb[offset-yoff+3+1]
								-rgb[offset-yoff-3]+1)));

						mlx += fabsf(0.1140*(rgb[offset+3+2]
								-rgb[offset-3+2]
								
								+0.5*(rgb[offset+yoff+3+2]
								-rgb[offset+yoff-3+2])
								
								+0.5*(rgb[offset-yoff+3+2]
								-rgb[offset-yoff-3+2])));


						float mly = fabsf(0.2989*(rgb[offset+yoff]
								-rgb[offset-yoff]

								+0.5*(rgb[offset+1+yoff]
								-rgb[offset+1-yoff])

								+0.5*(2*rgb[offset-1+yoff]
								-rgb[offset-1-yoff])
								));

						mly += fabsf(0.5870*(rgb[offset+yoff+1]
								-rgb[offset-yoff+1]

								+0.5*(rgb[offset+1+yoff+1]
								-rgb[offset+1-yoff+1])

								+0.5*(rgb[offset-1+yoff+1]
								-rgb[offset-1-yoff+1])
								));

						mly += fabsf(0.5870*(rgb[offset+yoff+2]
								-rgb[offset-yoff+2]

								+0.5*(rgb[offset+1+yoff+2]
								-rgb[offset+1-yoff+2])

								+0.5*(rgb[offset-1+yoff+2]
								-rgb[offset-1-yoff+2])
								));
						
							
						float k = (mlx+mly);
						f+=k;
					}
				}
				F[j]=f/averageI;//divide by overall intensity to get actual shape measure
				printf("%f-dist%f-",F[j],filmDist);
				delete [] rgb;
				sensor.WriteImage(1.f);

				delete[] samples;
				delete[] rays;
				delete[] Ls;
				delete[] Ts;
				delete[] isects;
			}

			int maxi=1;
			float maxf=0;
			//find last local maximum
			for(int k = 1;k < slice; k++){
				if(F[k]>F[k-1] && F[k]>F[k+1] && F[k]>maxf){
					maxi = k;
					maxf= F[k];
				}
			}
			UpdateDistance(filmDistStart+maxi*increment);
			filmDistEnd = filmDistStart +(maxi+1)*increment;
			filmDistStart = filmDistStart + (maxi-1)*increment;
			filmdists[i]=filmDist;
		}
		slice/=2;
	}
	filmDist= filmDistStart;
	int size = afZones.size();

	for(int m = 0; m < size; m++){
		if(filmdists[m]>filmDist)UpdateDistance(filmdists[m]);
	}
}