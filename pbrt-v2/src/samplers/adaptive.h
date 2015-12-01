
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SAMPLERS_ADAPTIVE_H
#define PBRT_SAMPLERS_ADAPTIVE_H

// samplers/adaptive.h*
#include "pbrt.h"
#include "sampler.h"
#include "../libsvm/svm.h"

//This means we actually initialize to 4+4 samples
#define ML_MIN_SAMPLES 4
//This gives us max 1024 samples: 4*2^8=1024
#define ML_MAX_LAYERS 2
//used in training
#define PIXEL_DIFF_THRESHOLD 0.01
const std::string DATA_PATH = "./samplers/mldata/";

typedef struct{
    float X,Y,Z;
    float var_Y, mean_Y;
    float delta_Y;//first half and second half
    bool single_shape;
    bool single_primitive;
    bool single_light;
    float var_ndepth;
    float *hull;//specified as for qhull
    int n_samples;
}pixel_data;

// AdaptiveSampler Declarations
class AdaptiveSampler : public Sampler {
public:
    // AdaptiveSampler Public Methods
    AdaptiveSampler(int xstart, int xend, int ystart, int yend, int ylength, int id_offset,
        int minSamples, int maxSamples, const string &type,
        float sopen, float sclose);
    Sampler *GetSubSampler(int num, int count);
    ~AdaptiveSampler();
    int RoundSize(int size) const {
        return RoundUpPow2(size);
    }
    int MaximumSampleCount() { return maxSamples; }
    int GetMoreSamples(Sample *sample, RNG &rng);
    bool ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
private:
    // AdaptiveSampler Private Methods
    bool needsSupersampling(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
    void combinePixelData();
    //Called for every training scene
    void writePixelToFile();
    void writeToFile(string filename, float value);
    //Called when triggered from script of last raining scene
    inline int get_id(){return idOffset + yLength*xPos + yPos;}
private:
    // AdaptiveSampler Private Data
    int xPos, yPos;
    int id;
    int yLength;
    int minSamples, maxSamples, currSamples;
    int idOffset;
    float *sampleBuf;
    enum MlType {  TRAIN, TEST};
    MlType datatype;
    int svmLayer;//-1 means stop
    pixel_data currPixel;
    pixel_data buffPixel;

};

void LoadModel();
int GetAndUpdateIDOffset(int newdata);
AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
    const Camera *camera);

#endif // PBRT_SAMPLERS_ADAPTIVE_H
