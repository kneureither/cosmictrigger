//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_PATTERNENGINE_H
#define COSMICTRIGGER_PATTERNENGINE_H

#include <vector>
#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <assert.h>
#include "SuperPixelHit.h"


class PatternEngine {

private:
    int getLayer(float r);
    int getXSP(float x);
    int getZSP(float z);
    int getPhiSP(float phi);

    template <typename T>
    int binSearch(std::vector<T> targetVector, T value, int start, int end);

    float getRadius(float x, float y);
    float getHitPhi(float x, float y);

    void initializeSPbins();

    float layerFactor[4] = {0.0}; // scale down x boundaries with reference to outer layer
    const float layerBoundaries[5] = {0.0, 26.00, 51.0, 78.5, 100.00}; //radial decision boundaries for layers

    int mode = 0; //default is cylindrical mode
    float wBinCount; //width (phi or x)
    float zBinCount;

    //only center region implemented for now
    float centralDetectorZmin = -200.00;
    float centralDetectorZmax = 200.00;
    float centralDetectorLength;

    float centralDetectorXmin; //etc

    std::vector<float> phiBins;
    std::vector<float> zBins;
    std::vector<float> xBins;

    std::string plottingpath = "../plots/PatternEngine/";

public:
    PatternEngine(float spXpartition, float spZpartition, int mode);
    PatternEngine(float spXpartition, float spZpartition, int mode, std::string plottingpath);
    unsigned int getSuperPixel(float x, float y, float z);
    int translateToSensorRefSP(int SPID);
    int translateToIndexRefSP(int SPID);
    void displayBinBoundaries();
    void testbinSearch();
};
#endif //COSMICTRIGGER_PATTERNENGINE_H
