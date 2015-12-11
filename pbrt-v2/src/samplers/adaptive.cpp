
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
                     float sopen, float sclose,svm_model **models,std::string modelpath)
    : Sampler(xstart, xend, ystart, yend, RoundUpPow2(max(mins, maxs)),
              sopen, sclose) {
    if(models == NULL)LoadModels(modelpath);
    else svmModels = models;
    xPos = xPixelStart;
    yPos = yPixelStart;
    yLength = ylength;
    idOffset = id_offset;
    svmModels = models;
    for(int f = 0; f <9; f++){
        buffPixel.nodes[f].index = f+1;
    }
    buffPixel.nodes[END].index = -1;

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
        //if(maxs > maxfromSVM)maxSamples = maxfromSVM;
        //else if(maxSamples < minSamples)maxSamples = minSamples;
        //else if (!IsPowerOf2(maxs)) maxSamples = RoundUpPow2(maxs)/2;
        //else maxSamples = maxs;
        maxSamples = maxfromSVM;
        //std::cout <<minSamples<<"*"<<maxSamples<<std::endl;
    }

    currSamples = ML_MIN_SAMPLES/2;
    sampleBuf = NULL;
}


AdaptiveSampler::~AdaptiveSampler() {
    delete[] sampleBuf;
    delete[] XYZBuf;
}


Sampler *AdaptiveSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new AdaptiveSampler(x0, x1, y0, y1, yLength, idOffset, minSamples, maxSamples, datatype == TRAIN ? "train" : "test", shutterOpen, shutterClose,svmModels,"");
}

//REWRITE
int AdaptiveSampler::GetMoreSamples(Sample *samples, RNG &rng) {
    if (!sampleBuf){
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples,maxSamples)];
        XYZBuf = new float[3*maxSamples];
    }
    if (svmLayer > 0) {
    //svmlayer = layer it has to passthrough AFTER this sample increase. i.e. layer 1 means it already has 2 rounds of samples of size min/2.
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
    
    if(datatype == TEST && currSamples == maxSamples){
        float layer = ML_MAX_LAYERS;
        writeToFile(DATA_PATH+"prediction",&layer,1);
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
        XYZBuf[3*i] *= XYZ_SCALING;
        XYZBuf[3*i+1] *= XYZ_SCALING;
        XYZBuf[3*i+2] *= XYZ_SCALING;
    }
    //init buff pixel
    buffPixel.nodes[DELTA_X].value = 0;
    buffPixel.nodes[DELTA_Y].value = 0;
    buffPixel.nodes[DELTA_Z].value = 0;
    buffPixel.nodes[M_X].value = 0;
    buffPixel.nodes[M_Y].value = 0;
    buffPixel.nodes[M_Z].value = 0;
    buffPixel.nodes[VAR_Y].value = 0;
    //buffPixel.nodes[DEPTH] = 0;
    buffPixel.nodes[SINGLESHAPE].value = true;

    float mX = XYZBuf[3*rangemin];
    float mY = XYZBuf[3*rangemin+1];
    float mZ = XYZBuf[3*rangemin+2];
    //float mR = RGBBuf[3*rangemin];
    //float mG = RGBBuf[3*rangemin+1];
    //float mB = RGBBuf[3*rangemin+2];
    float varY = 0;

    //calculate/analyze things
    for (int i = rangemin; i < rangemax -1; ++i){
        //shapes and primitive
        if (isects[i].shapeId != isects[i+1].shapeId)buffPixel.nodes[SINGLESHAPE].value = false;
        //color mean
        mX += XYZBuf[3*(i+1)];
        mY += XYZBuf[3*(i+1)+1];
        mZ += XYZBuf[3*(i+1)+2];
    }
    mX /= rangesize;
    mY /= rangesize;
    mZ /= rangesize;
    buffPixel.nodes[M_X].value = mX;
    buffPixel.nodes[M_Y].value = mY;
    buffPixel.nodes[M_Z].value = mZ;

    //try recursive once code is finished
    for (int i = rangemin; i < rangemax; ++i){
            varY += powf(XYZBuf[3*(i+1)+1] - mY,2);
    }

    varY /= rangesize;
    buffPixel.nodes[VAR_Y].value = varY;
    combinePixelData();

    if(svmLayer == 0)return true;
    if(datatype == TRAIN){
        writePixelToFile();
        if(svmLayer < ML_MAX_LAYERS)return true;
        else return false;
    }
    if(currSamples < minSamples)return true;
    double classification = svm_predict(svmModels[svmLayer-1], currPixel.nodes);
    float layer = svmLayer;
    if(classification <= 0) writeToFile(DATA_PATH+"prediction",&layer,1);
    return classification > 0;
}

void AdaptiveSampler::combinePixelData(){
    if(svmLayer == 0){
        //init currPixel with buffPixel
        currPixel = buffPixel;
    }else{
        currPixel.nodes[SINGLESHAPE].value = (currPixel.nodes[SINGLESHAPE].value && buffPixel.nodes[SINGLESHAPE].value);
        currPixel.nodes[DELTA_X].value = fabsf(currPixel.nodes[M_X].value - buffPixel.nodes[M_X].value);
        currPixel.nodes[DELTA_Y].value = fabsf(currPixel.nodes[M_Y].value - buffPixel.nodes[M_Y].value);
        currPixel.nodes[DELTA_Z].value = fabsf(currPixel.nodes[M_Z].value - buffPixel.nodes[M_Z].value);
        currPixel.nodes[VAR_Y].value = 0.5*(currPixel.nodes[VAR_Y].value + buffPixel.nodes[VAR_Y].value + 0.5*powf(currPixel.nodes[DELTA_Y].value,2));
        currPixel.nodes[M_X].value = 0.5*(currPixel.nodes[M_X].value + buffPixel.nodes[M_X].value);
        currPixel.nodes[M_Y].value = 0.5*(currPixel.nodes[M_Y].value + buffPixel.nodes[M_Y].value);
        currPixel.nodes[M_Z].value = 0.5*(currPixel.nodes[M_Z].value + buffPixel.nodes[M_Z].value);
    }
}

void AdaptiveSampler::writePixelToFile(){
    //debugging and visualizing
    writeToFile("DeltaX" + std::to_string(svmLayer), &(currPixel.nodes[DELTA_X]),1);
    writeToFile("DeltaY" + std::to_string(svmLayer), &(currPixel.nodes[DELTA_Y]),1);
    writeToFile("DeltaZ" + std::to_string(svmLayer), &(currPixel.nodes[DELTA_Z]),1);
    writeToFile("VarY" + std::to_string(svmLayer), &(currPixel.nodes[VAR_Y]),1);
    writeToFile("SingleShape" + std::to_string(svmLayer), &(currPixel.nodes[SINGLESHAPE]),1);
    //actual use
    float xyz[3] = {currPixel.nodes[M_X].value, currPixel.nodes[M_Y].value, currPixel.nodes[M_Z].value};
    writeToFile("PixelXYZ" + std::to_string(svmLayer), &xyz[0],3);
    writeToFile("AllRawData"+ std::to_string(svmLayer), currPixel.nodes,8);
}

void AdaptiveSampler::writeToFile(string filename, float* v, int size){
    std::ofstream file(DATA_PATH + filename, std::ios::app);
    if(file.is_open ()){
        file << id ; 
        for(int i = 0 ; i < size; i++){
            file << " "<< v[i];
        }
        file << "\n";
    }else Warning("Could not write to \"%s\"", filename.c_str());
    file.close();
}

void AdaptiveSampler::writeToFile(string filename, svm_node* v, int size){
    std::ofstream file(DATA_PATH + filename, std::ios::app);
    if(file.is_open ()){
        file << id ; 
        for(int i = 0 ; i < size; i++){
            file << " "<< v[i].value;
        }
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
        if(update < datacount)Warning("Might have overflow in nodes count");
        write << update;
    }else Warning("Could not write to count file");
    write.close();
    return datacount;
}

void AdaptiveSampler::LoadModels(std::string modelpath){
    svmModels = new svm_model*[ML_MAX_LAYERS];
    for(int i = 0; i < ML_MAX_LAYERS-1;i++){
        std::string s = modelpath+std::to_string(i+1);
        svmModels[i] = svm_load_model(s.c_str());
    }
}

AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int minsamp = params.FindOneInt("minsamples", ML_MIN_SAMPLES);
    int maxsamp = params.FindOneInt("maxsamples", 1024);
    if (PbrtOptions.quickRender) { minsamp = 2; maxsamp = 4; }
    string type = params.FindOneString("datatype", "train");
    string modelpath = params.FindOneString("modelpath", "");
    int ylength = yend - ystart;
    //need to deal with xstart, ystart !=0
    int id_offset = GetAndUpdateIDOffset(xend*yend);
    return new AdaptiveSampler(xstart, xend, ystart, yend, ylength, id_offset, minsamp, maxsamp, type,
         camera->shutterOpen, camera->shutterClose,NULL, modelpath);
}
