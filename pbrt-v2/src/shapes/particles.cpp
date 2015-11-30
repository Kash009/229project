/*#include "rng.h"
#include "particles.h"
#include "shape.h"
#include "shapes/sphere.h"
#include "geometry.h"
#include "paramset.h"
#include "accelerators/bvh.h"
#include "core/intersection.h"
#include <math.h>

  Particles::Particles(const Transform *o2w, const Transform *w2o, bool ro, 
        Point center, Point offset, float radius, int n_particles):
  Shape(o2w, w2o, ro), 
  		o2w(o2w),w2o(w2o),ro(ro),center(center),
  		offset(offset),radius(radius),n_particles(n_particles){};
 Particles::~Particles(){
    delete[] psamples;
    }
  BBox Particles::ObjectBound() const{
  	return BBox(Point(-radius,0,-radius),Point(-radius, 50, -radius));
  }
  BBox Particles::WorldBound() const{
  	return BBox(Point(center+Vector(-radius,0,-radius)),
  		Point(center +Vector(radius,50, radius)));
  }
    
   void Particles::Refine(vector<Reference<Shape> > &refined) const{
   		if(refined.size()==0)return;
   		UnsafeSampling(refined);

   		for (int i = 0; i < n_particles; ++i)
   		{
    
   			refined.push_back(
          new Sphere(new Transform(Translate(Vector(psamples[i])))*o2w),
          new Transform(Translate(-1*Vector(psamples[i])))*w2o),ro,0.1,-0.1,0.1,360));
   		}
   }

void Particles::UnsafeSampling(vector<Reference<Shape> > &refined)const{
   		if(refined.size()==0)return;
      BVHAccel *accel = new BVHAccel(&refined);
   		RNG rng = RNG();
   		psamples = new Point[n_particles];
   		int i = 0;
   		int j = 0;
   		while(i < n_particles && j <2*n_particles){
   			float r = rng.RandomFloat();
   			float theta = rng.RandomFloat();
   			float x = ((1-r)*offset.x + r*center.x) + radius*r*cosf(2*M_PI*theta);
   			float y = ((1-r)*offset.y + r*center.y) + radius*r*sinf(2*M_PI*theta);
   			Ray ray(Point(x,50,y), Vector(0,-1,0),0);
			Intersection isect;
			if(accel.Intersect(ray,&isect)){
				psamples[i] = ray(isect.rayEpsilon);
				i++;
			}
			j++;
		}
      delete accel;
   }

Particles *GenerateParticles(const Transform *o2w, const Transform *w2o, 
    bool ReverseOrientation, const ParamSet &params){
	int check;
    const float *c = params.FindFloat("center", &check);
    const float *o = params.FindFloat("offset", &check);
    float radius = params.FindOneFloat("radius", 1.f);
    int np = params.FindOneInt("particles", 27);
    return new Particles(o2w, w2o, ReverseOrientation, 
        Point(c[0],c[1],c[2]), Point(o[0],o[1],o[2]), radius, np);
}*/