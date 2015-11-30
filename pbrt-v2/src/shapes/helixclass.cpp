#include "stdafx.h"
#include "helixclass.h"
#include "shape.h"
#include "shapes/trianglemesh.h"
#include "geometry.h"
#include "paramset.h"
#include <math.h>

Helix::Helix(const Transform *o2w, const Transform *w2o, bool ro, 
        float thickness, unsigned int ridge, int n_twirls, 
        float twist, float r_start, bool orientation): 
		Shape(o2w, w2o, ro), 
		o2w(o2w), w2o(w2o),ro(ro),
		thickness(thickness), xsect_nridge(ridge), 
		twirls(n_twirls), twist(twist), r_start(r_start),twirl_orient(orientation){
		sample_rate = 256;
		generateSkeleton();
}

Helix::~Helix(){
	delete[] skeleton; 
} 

float Helix::d(float t){
	return r_start*(1 - powf(t/(float)twirls,3));
}

float crsct(float t, int ridge){
	return 0.76+0.12*(powf((sin(ridge*M_PI*t)),4)+powf((sin(ridge*M_PI*t)),2))
	+0.015*(powf((sin(4.7*ridge*M_PI*t)),4)+powf((sin(4.7*ridge*M_PI*t)),2));
}

void Helix::generateSkeleton(){
	int orient = 1;
	if(twirl_orient)orient = -1;

	skeleton = new HelixSample[sample_rate*twirls];
	float increment = 1.0f/sample_rate;
	float wrinkle = 1/(float)twirls;	
	//first point
	skeleton[0].center = Point(d(increment), 0, 0);
	skeleton[0].radius = 0;
	skeleton[0].dir = Vector(0,orient*skeleton[0].center.x,0);
	skeleton[0].up = Vector(0,0,1);
	skeleton[0].theta= 0;
	//second point
	skeleton[1].center = Point(d(increment)*cosf(orient*2*M_PI*increment), d(increment)*sinf(orient*2*M_PI*increment), 0.8*increment*thickness);
	skeleton[1].radius = increment*thickness;
	skeleton[1].theta= increment*twist;
	
	//rest of first twirl, with direction defined by vector connecting previous and next points
	for(int i = 2; i < sample_rate; i++){
		skeleton[i].center=Point(orient*d(i*increment)*cosf(2*M_PI*i*increment),
							 orient*d(i*increment)*sinf(2*M_PI*i*increment),
							 0.8*i*increment*thickness);
		skeleton[i].radius= i*increment*thickness;
		skeleton[i].theta= i*increment*twist;
		Vector dir = Normalize(Vector(skeleton[i].center)- Vector(skeleton[i-2].center));
		skeleton[i-1].dir = dir;
		//Vector whose plane spanned by itself and dir includes (0.0.1),
			//i.e. the vector equivalent to "up" in this space
		skeleton[i-1].up = Normalize(Vector(-dir.x*dir.z,-dir.y*dir.z, dir.x*dir.x + dir.y*dir.y));
	}

	//twirl dependent on first
	if(twirls < 2)return;
	else if(twirls >2){
		for(int i = 0; i < sample_rate; i++){
			int k = sample_rate + i;
			float h = powf((1+i*increment)*thickness,2) - powf(d(k*increment)-d((k-sample_rate)*increment),2);
			if(h<0)h=0;
			skeleton[k].center=Point(orient*d(k*increment)*cosf(2*M_PI*i*increment), 
									orient*d(k*increment)*sinf(2*M_PI*i*increment),
									skeleton[k-sample_rate].center.z+(1.0+0.04*(1-i*increment))*0.8*sqrtf(h));
			skeleton[k].radius= thickness*(0.95+(1-k*increment/(float)(twirls-1))
								*(0.01*cosf(100*M_PI*i*increment)
								+0.02*cosf(0.5+78*M_PI*i*increment+1)
								+0.06*cosf(1+22*M_PI*i*increment+1)));
			skeleton[k].theta= k*increment*twist;
			Vector dir = Normalize(Vector(skeleton[k].center)- Vector(skeleton[k-2].center));
			skeleton[k-1].dir = dir;
			skeleton[k-1].up = Normalize(Vector(-dir.x*dir.z,-dir.y*dir.z, dir.x*dir.x + dir.y*dir.y));
		}
	}
	//middle twirls
	for(int j = 2; j < twirls-1;j++){
		for(int i = 0; i < sample_rate; i++){
			int k = j*sample_rate + i;
			float h = powf(2*thickness,2) - powf(d(k*increment)-d((k-sample_rate)*increment),2);
			if(h<0)h=0;
			skeleton[k].center=Point(orient*d(k*increment)*cosf(2*M_PI*i*increment), 
									orient*d(k*increment)*sinf(2*M_PI*i*increment),
									skeleton[k-sample_rate].center.z+0.8*sqrtf(h));
			if(h < 0.001)skeleton[k].center.z = 2.5*skeleton[k-1].center.z-1.5*skeleton[k-2].center.z;
			skeleton[k].radius= thickness*(0.95+(1-k*increment/(float)(twirls-1))*
								(0.01*cosf(130*(1-j*wrinkle)*M_PI*i*increment)
								+0.02*cosf(0.5+78*(1-j*wrinkle)*M_PI*i*increment+1))
								+0.07*cosf(1+44*(1-j*wrinkle)*M_PI*i*increment+1));
			skeleton[k].theta= k*increment*twist;
			Vector dir = Normalize(0.1*skeleton[k-2].dir
								+0.8*Normalize(Vector(skeleton[k].center)- Vector(skeleton[k-2].center)));
			skeleton[k-1].dir = dir;
			skeleton[k-1].up = Normalize(Vector(-dir.x*dir.z,-dir.y*dir.z, dir.x*dir.x + dir.y*dir.y));
		}
	}
	float init = skeleton[(twirls-1)*sample_rate-2].center.z - skeleton[(twirls-1)*sample_rate-3].center.z;
	//last swirl, twirls upwards fairly sharply
	for(int x = 0; x < sample_rate; x++){
		int y = (twirls-1)*sample_rate + x;
		float h;
		if(x==0)h = 1.82*thickness;
		else h = 1.82*thickness*(1+2*powf(x*increment,7));
		if(h<0)h=0;
		skeleton[y].center=Point(orient*d(y*increment)*cosf(2*M_PI*x*increment), 
								orient*d(y*increment)*sinf(2*M_PI*x*increment),
								skeleton[y-sample_rate].center.z + 0.8*h);	
		skeleton[y].radius = thickness*(1-powf(x*increment,8));//*(0.98+0.02*cosf(20*M_PI*x*increment));
		skeleton[y].theta= twist*increment*(y-powf(increment*x,0.5));
		Vector dir = Normalize((powf(increment*x,0.5))*Vector(0,0,1)+(1-powf(increment*x,0.5))*skeleton[y-sample_rate].dir);
		skeleton[y-1].dir = Normalize(dir);
		skeleton[y-1].up = Normalize(Vector(-dir.x*dir.z,-dir.y*dir.z, dir.x*dir.x + dir.y*dir.y));
	}
}

void Helix::Refine(vector<Reference<Shape> > &refined) const{
	int xsect_sample = 256;
	//int xsect_nridge = 6;
	int per_xsect = xsect_sample * xsect_nridge;
	float *uv = new float[2*(per_xsect*(sample_rate*twirls -2))];
	Point *P = new Point[per_xsect*(sample_rate*twirls -2) + 2];
	int *V = new int[3*per_xsect*(sample_rate*twirls -2)];
	Normal *N= new Normal[per_xsect*(sample_rate*twirls -2) + 2];
	
	//UVs
	uv[0] = 0.5;
	uv[1] = 0.0;
	for(int i = 0;i < sample_rate*twirls -2 ; i++){
		for(int j = 1; j <= per_xsect; j++){
					uv[2*(per_xsect*i+j)]=	40*j/(float)per_xsect;
					uv[2*(per_xsect*i+j)+1]=4*(i%sample_rate)/(float)(sample_rate);
		}
	}
	uv[2*(per_xsect*(sample_rate*twirls -2) + 1)]= 0.5;
	uv[2*(per_xsect*(sample_rate*twirls -2) + 1)+1] = 1.0;
	//P
	float rad[per_xsect];
	for(int i = 0; i < per_xsect; i++){
			rad[i] = crsct(i/(float)per_xsect,xsect_nridge);
	}
	int index = 0;
	P[index++]=skeleton[0].center;
	for(int i = 0;i < sample_rate*twirls -2 ; i++){
		Vector x = Normalize(Cross(skeleton[i].up,skeleton[i].dir));
		Vector y = Normalize(skeleton[i].up);
		Point p = skeleton[i].center;
		float t = -1*skeleton[i].theta;
		float r = skeleton[i].radius;
		for(int j = 0; j < per_xsect; j++){
			P[index ++]=p + r*rad[j]*(cosf(2*M_PI*(t+j/(float)per_xsect))*x 
								+ sinf(2*M_PI*(t+j/(float)per_xsect))*y);
		}
	}
	P[index++]=skeleton[sample_rate*twirls -1].center;
	//N
	N[0]=Normalize(Normal(-1*skeleton[0].dir));
	for(int i = 0;i < sample_rate*twirls -2 ; i++){
		N[i*per_xsect]= Normalize(Normal(2*Vector(P[i*per_xsect])
								-Vector(P[i*per_xsect +1])
								-Vector(P[i*per_xsect + per_xsect])));
		for(int j = 2; j <= per_xsect; j++){
			N[i*per_xsect+j]= Normalize(Normal(2*Vector(P[i*per_xsect])
									-Vector(P[i*per_xsect + (j+1)%per_xsect])
									-Vector(P[i*per_xsect + (j-1)%per_xsect])));
		}
	}
	N[per_xsect*(sample_rate*twirls-2)+1]= Normalize(Normal(skeleton[sample_rate*twirls-1].dir));
	//V
	index = 0;
	for(int i = 1; i <=per_xsect;i++){
		V[index++]=0;
		V[index++]=i;
		V[index++]=1+(i%per_xsect);
	}
	for(int i = 0; i< twirls*sample_rate-3;i++){
		for(int j = 1; j <=per_xsect;j++){
			V[index++]=per_xsect*i + j;
			V[index++]=per_xsect*(i+1) + j;
			V[index++]=per_xsect*(i+1) + 1+(j%per_xsect);

			V[index++]=per_xsect*i + 1+(j%per_xsect);
			V[index++]=per_xsect*i + j;
			V[index++]=per_xsect*(i+1) +1+(j%per_xsect);
		}
	}

	for(int i = 1; i <=per_xsect;i++){
		int j = (per_xsect*(sample_rate*twirls-3));
		V[index++]=1+j+(i%per_xsect);
		V[index++]=j+i;
		V[index++]=per_xsect*(sample_rate*twirls -2)+1;
	}

	refined.push_back(new TriangleMesh(o2w, w2o, ro, 
						2*per_xsect*(sample_rate*twirls-2), 
						2+per_xsect*(sample_rate*twirls-2),
						 V, P, NULL, NULL, uv, NULL));
	delete[] P;
	delete[] V;
	delete[] uv;
	delete[] N;

}

BBox Helix::ObjectBound() const {
    return BBox(Point(-(thickness+r_start), -(thickness+r_start), 0),
                Point( (thickness+r_start),(thickness+r_start), thickness + skeleton[sample_rate*twirls-1].center.z));
}

BBox Helix::WorldBound() const{
	Point p1 = (*ObjectToWorld)(Point(-(thickness+r_start), -(thickness+r_start), 0));
	Point p2 = (*ObjectToWorld)(Point( (thickness+r_start),(thickness+r_start), thickness + skeleton[sample_rate*twirls-1].center.z));
    BBox box = Union(box, p1);
    box = Union(box, p2);
    return box;
}

Helix *CreateHelixShape(const Transform *o2w, const Transform *w2o, 
    bool ReverseOrientation, const ParamSet &params){
    
    float base = params.FindOneFloat("radius", 1.f);
    int n_twirls = params.FindOneInt("twirls", 4);
	float thickness = params.FindOneFloat("thickness", 0.5f*base);
	int ridge = params.FindOneInt("ridge", 6);
	float twist = params.FindOneFloat("twist", 1.0f);
	bool orientation = !(params.FindOneInt("orientation",0));
    return new Helix(o2w, w2o, ReverseOrientation, thickness,
                      ridge, n_twirls, twist, base, orientation);
}
