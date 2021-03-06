
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
    nodes, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
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
#include "svm.h"

//This means we actually initialize to 4+4 samples
#define ML_MIN_SAMPLES 8
//This gives us max 1024 samples: 8*2^7=512
#define ML_MAX_LAYERS 6
//scale values back
#define XYZ_SCALING 100
const std::string DATA_PATH = "../data/";

enum Feature{ DELTA_X, DELTA_Y, DELTA_Z, M_X, M_Y, M_Z, VAR_Y,/*DEPTH,*/ SINGLESHAPE,END};

typedef struct{
    svm_node nodes[9];
}pixel_data;

// AdaptiveSampler Declarations
class AdaptiveSampler : public Sampler {
public:
    // AdaptiveSampler Public Methods
    AdaptiveSampler(int xstart, int xend, int ystart, int yend, int ylength, int id_offset,
        int minSamples, int maxSamples, const string &type,
        float sopen, float sclose, svm_model** models,std::string modelpath);
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
    void LoadModels(std::string modelpath);
    bool needsSupersampling(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
    void combinePixelData();
    //Called for every training scene
    void writePixelToFile();
    void writeColor(string filename, float x, float y, float z);
    void writeToFile(string filename, float *value, int datasize);
    void writeToFile(string filename, svm_node* value, int datasize);
    //Called when triggered from script of last raining scene
    inline int get_id(){return idOffset + yLength*xPos + yPos;}
private:
    // AdaptiveSampler Private nodes
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
    svm_model** svmModels;
};

int GetAndUpdateIDOffset(int newdata);
AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
    const Camera *camera);

#endif // PBRT_SAMPLERS_ADAPTIVE_H
