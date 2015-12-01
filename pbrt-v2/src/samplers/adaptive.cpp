
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


// samplers/adaptive.cpp*
#include "stdafx.h"
#include "samplers/adaptive.h"
#include "paramset.h"
#include "film.h"
#include "primitive.h"
#include "intersection.h"
#include "camera.h"
#include "montecarlo.h"
#include <iostream>
#include <fstream>

//Call method
//while(getmoresamples){if(reportresults)writepixeltofile}
// AdaptiveSampler Method Definitions
AdaptiveSampler::AdaptiveSampler(int xstart, int xend,
                     int ystart, int yend, int ylength, int id_offset, int mins, int maxs, const string &t,
                     float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, RoundUpPow2(max(mins, maxs)),
              sopen, sclose) {
    
    xPos = xPixelStart;
    yPos = yPixelStart;
    yLength = ylength;
    idOffset = id_offset;

    if(t == "train") datatype = TRAIN;
    else datatype = TEST;
    //check for redundancy when more awake
    minSamples = ML_MIN_SAMPLES;
    currSamples = minSamples/2;
    svmLayer = 0;

    if(datatype == TEST){
        if (maxs < ML_MIN_SAMPLES){
            Warning("Maximum samples has to be greater or equal to %d; setting to %d", 
                ML_MIN_SAMPLES, ML_MIN_SAMPLES);
            maxSamples = ML_MIN_SAMPLES;

        }else if (!IsPowerOf2(maxs)) {
            Warning("Maximum pixel samples being rounded down to power of 2");
            maxSamples = RoundUpPow2(maxs)/2;
            //remember to truncate to ml max layers?
        }else
            maxSamples = maxs;
    }else{
        //for training
        //doublecheck redundancy
        maxSamples = ML_MIN_SAMPLES * powf(2,ML_MAX_LAYERS);
    }
    sampleBuf = NULL;
}


AdaptiveSampler::~AdaptiveSampler() {
    delete[] sampleBuf;
}


Sampler *AdaptiveSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new AdaptiveSampler(x0, x1, y0, y1, yLength, idOffset, minSamples, maxSamples, datatype == TRAIN ? "train" : "test", shutterOpen, shutterClose);
}

//REWRITE
int AdaptiveSampler::GetMoreSamples(Sample *samples, RNG &rng) {
    if (!sampleBuf)
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples,maxSamples)];
    if (svmLayer > 0) {
    //svmlayer = layer it has to passthrough AFTER this sample increase. i.e. layer 1 means it already has 2 rounds of samples.
        LDPixelSample(xPos, yPos, shutterOpen, shutterClose, currSamples,
                      &(samples[currSamples]), sampleBuf, rng);
        currSamples *= 2;
        return currSamples;
    }
    //end of image
    else if(svmLayer < 0){
        if(yPos == yPixelEnd) return 0;
        svmLayer = 0;//set to new pixel
    }
    //new pixel
    LDPixelSample(xPos, yPos, shutterOpen, shutterClose, currSamples, samples, sampleBuf, rng);
    //note that we're not incrementing currsamples for the base case
    return currSamples;
}


bool AdaptiveSampler::ReportResults(Sample *samples,
        const RayDifferential *rays, const Spectrum *Ls,
        const Intersection *isects, int count) {
  
    if(needsSupersampling(samples, rays, Ls, isects, count)){
        svmLayer++;
        //increase number of samples
        return false;
    }
    if (++xPos == xPixelEnd) {
        xPos = xPixelStart;
        ++yPos;
    }
    svmLayer = -1;
    currSamples = minSamples/2;
    return true;
}


bool AdaptiveSampler::needsSupersampling(Sample *samples,
        const RayDifferential *rays, const Spectrum *Ls,
        const Intersection *isects, int count) {
    
    int rangemin;
    int rangemax;
    int rangesize;
    if(svmLayer == 0){
        id = get_id();
        rangemin = 0;
        rangemax = currSamples;
    }else{
        rangemin = currSamples/2;
        rangemax = currSamples;
    }
    rangesize = rangemax - rangemin;
    //parse things
    //init buff pixel
    buffPixel.single_shape = true;
    buffPixel.single_primitive = true;
    float mean = 0;
    //calculate/analyze things
    for (int i = rangemin; i < rangemax -1; ++i){
        if (isects[i].shapeId != isects[i+1].shapeId)buffPixel.single_shape = false;
        mean += Ls[i].y();
    }
    mean /= rangesize;
    combinePixelData();
    if(datatype == TRAIN){
        writePixelToFile();
        if(svmLayer < ML_MAX_LAYERS)return true;
    }
    //currently test doesn't do anything besides minimum samples
    if(svmLayer == 0)return true;
    return false;
}

void AdaptiveSampler::combinePixelData(){
    if(svmLayer == 0){
        //init currPixel with buffPixel
        currPixel = buffPixel;
    }else{
        currPixel.single_shape = (currPixel.single_shape && buffPixel.single_shape);
        currPixel.mean_Y = (currPixel.mean_Y + buffPixel.mean_Y)/2;
    }
}

void AdaptiveSampler::writePixelToFile(){
    writeToFile("SingleShape" + std::to_string(svmLayer), currPixel.single_shape);
}

void AdaptiveSampler::writeToFile(string filename, float value){
    std::ofstream file(DATA_PATH + filename, std::ios::app);
    if(file.is_open ()){
        file << id << " " << value << "\n";
    }else Warning("Could not write to \"%s\"", filename.c_str());
    file.close();
}

int GetAndUpdateIDOffset(int newpoints){
    std::ifstream read(DATA_PATH + "datacount");
    int datacount = 0;
    if(read.is_open()){
        read >> datacount;
    }else Warning("No pre-existing count file, setting to zero.");
    read.close();


    std::ofstream write(DATA_PATH + "datacount");
    if(write.is_open()){
        int update = datacount + newpoints;
        if(update < datacount)Warning("Might have overflow in data count");
        write << update;
    }else Warning("Could not write to count file");
    write.close();
    return datacount;
}

AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int minsamp = params.FindOneInt("minsamples", ML_MIN_SAMPLES);
    int maxsamp = params.FindOneInt("maxsamples", 64);
    if (PbrtOptions.quickRender) { minsamp = 2; maxsamp = 4; }
    string type = params.FindOneString("datatype", "train");
    int ylength = yend - ystart;
    //need to deal with xstart, ystart !=0
    int id_offset = GetAndUpdateIDOffset(xend*yend);
    return new AdaptiveSampler(xstart, xend, ystart, yend, ylength, id_offset, minsamp, maxsamp, type,
         camera->shutterOpen, camera->shutterClose);
}