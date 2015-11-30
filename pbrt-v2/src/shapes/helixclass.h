#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HELIX_H
#define PBRT_SHAPES_HELIX_H

#include "shape.h"

typedef struct{
    Point center;
    Vector dir;
    Vector up;
    float radius;
    float theta;
}HelixSample;

// Cylinder Declarations
class Helix : public Shape {
public:
    // Helix methods
    Helix(const Transform *o2w, const Transform *w2o, bool ro, 
        float thickness, unsigned int ridge, int n_twirls, 
        float twist, float r_start, bool orientation);
    ~Helix();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const {return false;}
    void Refine(vector<Reference<Shape> > &refined) const;

private:
    //Private Data
    const Transform *o2w;
    const Transform *w2o;
    bool ro;
    bool twirl_orient;
    HelixSample *skeleton;
    void generateSkeleton();
    float d(float t);
   // float crsct(float t);
    //float t_end;//is n_twirls + 1.5(a twirl for beginning, half for end)
    int xsect_nridge;
    float thickness;
    float r_start;
    int twirls;
    float twist;
    //float type;//for cream/softserve/etc
   
    //for sampling skeleton
    int sample_rate;
    //RandomSampler
};

Helix *CreateHelixShape(const Transform *o2w, const Transform *w2o, 
    bool ReverseOrientation, const ParamSet &params);
//Helix *printObj();
#endif

