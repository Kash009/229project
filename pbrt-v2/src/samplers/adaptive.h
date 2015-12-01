
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
#define ML_MIN_SAMPLES 8
//This gives us max 1024 samples: 8*2^7=1024
#define ML_MAX_LAYERS 6
//used in training // currently only difference in illuminance, can do more complicate things in Lab color space
#define PIXEL_DIFF_THRESHOLD 0.01
//values too small, offset rounding error
#define XYZ_SCALING 100
const std::string DATA_PATH = "./samplers/mldata/";

typedef struct{
    float X,Y,Z;//pixel value in ciexyz space
    float var_Y;//illuminance variance
    float mean_bounces;//n light Bounces
    float var_intersection;//TODO: determine quantification
    float delta_Y;//first half and second half
    bool single_shape;//object component
    bool single_primitive;//object
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
    void writeColor(string filename, float x, float y, float z);
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
    float *XYZBuf;
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
