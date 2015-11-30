/*#ifndef PBRT_SHAPES_HELIX_H
#define PBRT_SHAPES_HELIX_H

#include "shape.h"

// Cylinder Declarations
class Particles : public Shape {
public:
    // Helix methods
    Particles(const Transform *o2w, const Transform *w2o, bool ro, 
        Point center, Point offset, float radius, int n_particles);
    ~Particles();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const {return false;}
    void Refine(vector<Reference<Shape> > &refined) const;

private:
    //Private Data
    void UnsafeSampling(vector<Reference<Shape> > &refined) const;
    const Transform *o2w;
    const Transform *w2o;
    bool ro;
    const Point center;
    const Point offset;
    const float radius;
    const int n_particles;
    Point *psamples;
};

Particles *GenerateParticles(const Transform *o2w, const Transform *w2o, 
    bool ReverseOrientation, const ParamSet &params);

#endif*/

