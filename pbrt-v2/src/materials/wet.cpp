/*#include "stdafx.h"
#include "materials/wet.h"
#include "textures/constant.h"
#include "volume.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

BSDF *WetMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
              const DifferentialGeometry &dgShading,
              MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)  Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum R = Kr0->Evaluate(dgs).Clamp();
    float e = eta->Evaluate(dgs);
    if (!R.IsBlack()){
                BxDF* bxdf = BSDF_ALLOC(arena, FresnelBlend)(Kd0, Ks0, BSDF_ALLOC(arena, Blinn)());
                    bsdf->Add(bxdf);
                }
    return bsdf;
}


BSSRDF *WetMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom,
              const DifferentialGeometry &dgShading,
              MemoryArena &arena) const {

    float e = eta->Evaluate(dgShading);
    float mfp = meanfreepath->Evaluate(dgShading);
    Spectrum kd = Kd->Evaluate(dgShading).Clamp();
    Spectrum sigma_a, sigma_prime_s;
    SubsurfaceFromDiffuse(kd, mfp, e, &sigma_a, &sigma_prime_s);
    return BSDF_ALLOC(arena, BSSRDF)(sigma_a, sigma_prime_s, e);
}


WetMaterial *CreateWetMaterial(const Transform &xform,
        const TextureParams &mp) {
    float Kd[3] = { .5, .5, .5 };
    Reference<Texture<Spectrum> > kd0 = mp.GetSpectrumTexture("Kd0", Spectrum::FromRGB(Kd));
    Reference<Texture<Spectrum> > kd1 = mp.GetSpectrumTexture("Kd1", Spectrum::FromRGB(Kd));
    return new WetMaterial();
}
*/
