/*
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_WET_H
#define PBRT_MATERIALS_WET_H

// materials/kdsubsurface.h*
#include "pbrt.h"
#include "material.h"

// KdSubsurfaceMaterial Declarations
class WetMaterial : public Material {
public:
    // KdSubsurfaceMaterial Public Methods
    WetMaterial(Reference<Texture<Spectrum> > kd0,
            Reference<Texture<Spectrum> > kr0,
            float n0,
            Reference<Texture<Spectrum> > kd1,
            Reference<Texture<Spectrum> > kr1,
            float n1,
            Reference<Texture<float> > mfp,
            Reference<Texture<float> > e,
            Reference<Texture<float> > bump) {
        Kd0 = kd0;
        Kr0 = kr0;
        Kd1 = kd1;
        Kr1 = kr1;
        N0=n0;
        N1=n1;
        meanfreepath = mfp;
        eta = e;
        bumpMap = bump;
    }
    
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
private:
    Reference<Texture<Spectrum> > Kd0, Kr0, Kd1, Kr1;
    float N0,N1;
    Reference<Texture<float> > meanfreepath, eta, bumpMap;

};


WetMaterial *CreateWetMaterial(const Transform &xform,
        const TextureParams &mp);

#endif*/
