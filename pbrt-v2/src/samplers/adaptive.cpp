
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
#include <string>

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
    
    int maxfromSVM = ML_MIN_SAMPLES * powf(2,ML_MAX_LAYERS);
    
    svmLayer = 0;
    if(datatype == TRAIN){
        //for training
        //doublecheck redundancy
        minSamples = ML_MIN_SAMPLES;
        maxSamples = maxfromSVM;

    }else{
        minSamples = RoundUpPow2(mins);
        if(minSamples > maxfromSVM) minSamples = maxfromSVM;
        else if(minSamples < ML_MIN_SAMPLES) minSamples = ML_MIN_SAMPLES;
        if(maxs > maxfromSVM)maxSamples = maxfromSVM;
        else if(maxSamples < minSamples)maxSamples = minSamples;
        else if (!IsPowerOf2(maxs)) maxSamples = RoundUpPow2(maxs)/2;
        else maxSamples = maxs;
    }

    currSamples = ML_MIN_SAMPLES/2;
    sampleBuf = NULL;
}


AdaptiveSampler::~AdaptiveSampler() {
    delete[] sampleBuf;
    delete[] XYZBuf;
    delete[] RGBBuf;
}


Sampler *AdaptiveSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new AdaptiveSampler(x0, x1, y0, y1, yLength, idOffset, minSamples, maxSamples, datatype == TRAIN ? "train" : "test", shutterOpen, shutterClose);
}

//REWRITE
int AdaptiveSampler::GetMoreSamples(Sample *samples, RNG &rng) {
    if (!sampleBuf){
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples,maxSamples)];
        XYZBuf = new float[3*maxSamples];
        RGBBuf = new float[3*maxSamples];
    }
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
    currSamples = ML_MIN_SAMPLES/2;
    return true;
}


bool AdaptiveSampler::needsSupersampling(Sample *samples,
        const RayDifferential *rays, const Spectrum *Ls,
        const Intersection *isects, int count) {
    
    if(currSamples == maxSamples){
        if(datatype == TRAIN){

        }
        return false;
    }
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
    //
    //convert to XYZ. TODO: find documentation to stop initial xyz to rgb! maybe implement to Lab as well.
    for (int i = rangemin; i < rangemax; ++i){
        Ls[i].ToXYZ(&XYZBuf[3*i]);
        Ls[i].ToRGB(&RGBBuf[3*i]);
       // XYZBuf[3*i] *= XYZ_SCALING;
       // XYZBuf[3*i +1] *= XYZ_SCALING;
       // XYZBuf[3*i +2] *= XYZ_SCALING;
    }
    //init buff pixel
    buffPixel.data[DELTA_X] = 0;
    buffPixel.data[DELTA_Y] = 0;
    buffPixel.data[DELTA_Z] = 0;
    buffPixel.data[VAR_Y] = 0;
    //buffPixel.data[DEPTH] = 0;
    buffPixel.data[SINGLESHAPE] = true;

    float mBounces = rays[rangemin].depth;
    float mX = XYZBuf[3*rangemin];
    float mY = XYZBuf[3*rangemin+1];
    float mZ = XYZBuf[3*rangemin+2];
    float mR = RGBBuf[3*rangemin];
    float mG = RGBBuf[3*rangemin+1];
    float mB = RGBBuf[3*rangemin+2];
    float varY = 0;

    //calculate/analyze things
    for (int i = rangemin; i < rangemax -1; ++i){
        //shapes and primitive
        if (isects[i].shapeId != isects[i+1].shapeId)buffPixel.data[SINGLESHAPE] = false;
        //color mean
        mX += XYZBuf[3*(i+1)];
        mY += XYZBuf[3*(i+1)+1];
        mZ += XYZBuf[3*(i+1)+2];
        mBounces += rays[i+1].depth;
    }
    mX /= rangesize;
    mY /= rangesize;
    mZ /= rangesize;
    mR /= rangesize;
    mG /= rangesize;
    mB /= rangesize;
    float mrgb[3] = {mR,mG,mB};
    if(svmLayer == 0)writeToFile("rgb",mrgb,3);
    mBounces /=rangesize;

    buffPixel.XYZ[0] = mX;
    buffPixel.XYZ[1] = mY;
    buffPixel.XYZ[2] = mZ;

    //try recursive once code is finished
    for (int i = rangemin; i < rangemax; ++i){
            varY += powf(XYZBuf[3*(i+1)+1] - mY,2);
    }
    varY /= rangesize;
    buffPixel.data[VAR_Y] = varY;
    combinePixelData();

    if(svmLayer == 0)return true;
    if(datatype == TRAIN){
        writePixelToFile();
        if(svmLayer < ML_MAX_LAYERS)return true;
        else return false;
    }
    //writeToFile("sampleCount", currSamples);
    if(currSamples < minSamples)return true;

    return false;
}

void AdaptiveSampler::combinePixelData(){
    if(svmLayer == 0){
        //init currPixel with buffPixel
        currPixel = buffPixel;
    }else{
        currPixel.data[SINGLESHAPE] = (currPixel.data[SINGLESHAPE] && buffPixel.data[SINGLESHAPE]);
        currPixel.data[DELTA_X] = fabsf(currPixel.XYZ[0] - buffPixel.XYZ[0]);
        currPixel.data[DELTA_Y] = fabsf(currPixel.XYZ[1] - buffPixel.XYZ[1]);
        currPixel.data[DELTA_Z] = fabsf(currPixel.XYZ[2] - buffPixel.XYZ[2]);
        currPixel.data[VAR_Y] = 0.5*(currPixel.data[VAR_Y] + buffPixel.data[VAR_Y] + 0.5*powf(currPixel.data[DELTA_Y],2));
    
        currPixel.XYZ[0] = (currPixel.XYZ[0] + buffPixel.XYZ[0])/2;
        currPixel.XYZ[1] = (currPixel.XYZ[1] + buffPixel.XYZ[1])/2;
        currPixel.XYZ[2] = (currPixel.XYZ[2] + buffPixel.XYZ[2])/2;
    }
}

void AdaptiveSampler::writePixelToFile(){
    //debugging and visualizing
    writeToFile("DeltaX" + std::to_string(svmLayer), &(currPixel.data[DELTA_X]),1);
    writeToFile("DeltaY" + std::to_string(svmLayer), &(currPixel.data[DELTA_Y]),1);
    writeToFile("DeltaZ" + std::to_string(svmLayer), &(currPixel.data[DELTA_Z]),1);
    writeToFile("VarY" + std::to_string(svmLayer), &(currPixel.data[VAR_Y]),1);
    writeToFile("SingleShape" + std::to_string(svmLayer), &(currPixel.data[SINGLESHAPE]),1);
    //actual use
    writeToFile("AllRawData"+ std::to_string(svmLayer), currPixel.data,6);
    writeToFile("PixelXYZ" + std::to_string(svmLayer), currPixel.XYZ,3);
}

void AdaptiveSampler::writeToFile(string filename, float* value, int size){
    std::ofstream file(DATA_PATH + filename, std::ios::app);
    if(file.is_open ()){
        file << id ; 
        for(int i = 0 ; i < size; i++)
            file << " "<< value[i];
        file << "\n";
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
