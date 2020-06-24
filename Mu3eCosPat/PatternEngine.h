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


//TODO Implement other mapping of superpixels (maybe just a cartesian one
// -- This could either be done by deriving cartesian and phi dependent from a common base class
// -- Or include it into this class and use the int mode parameter


class PatternEngine {

private:
    int getLayer(const float r);
    int getXSP(const float x);
    int getZSP(const float z);
    int getPhiSP(float phi);
    int getSPZcoord(int sp2DID);
    int getSPXcoord(int sp2DID);

    template <typename T>
    int binSearch(std::vector<T> targetVector, T value, int start, int end);

    template <typename T>
    int seqSearch(std::vector<T> targetVector, T value);

    float getRadius(float x, float y);
    float getHitPhi(float x, float y);

    void initializeSPbins();
    void initializeMembers(float spXpartition, float spZpartition, int mode);

    std::string pltf();


    float layerFactor[4] = {0.0}; // scale down x boundaries with reference to outer layer
    const float layerBoundaries[5] = {0.0, 26.00, 51.0, 78.5, 100.00}; //radial decision boundaries for layers

    int mode = 0; //default is cylindrical mode
    int wBinCount; //width (phi or x)
    int zBinCount;
    int totalBinCount;

    //only center region implemented for now
    float centralDetectorZmin = -200.00;
    float centralDetectorZmax = 200.00;
    float centralDetectorLength;

    float widthMax;
    float widthMin;
    float widthDimLength;

    float centralDetectorXmin; //etc

    std::vector<float> phiBins;
    std::vector<float> zBins;
    std::vector<float> xBins;

    float SPZsize;
    float SPPhisize;

    std::vector<int> spWeights[4];

    std::string plottingpath = "../plots/PatternEngine/";
    std::string plottingfile = "PEplots.pdf";
    std::string runspecs;
    int firstplot = 1;

public:
    PatternEngine(float spXpartition, float spZpartition, int mode);
    PatternEngine(float spXpartition, float spZpartition, int mode, std::string plottingpath);
//    ~PatternEngine();

    unsigned int getSuperPixel(const float x, const float y, const float z);
    int getLayerFromSPID(const unsigned int SPID);
    int getAreaFromSPID(const unsigned int SPID);
    int getIndexFromSPID(const unsigned int SPID);
    int computeIndex(const int zSPIndex, const int phiSPIndex);
    unsigned int computeSPID(const int layer, const int area, const unsigned int index);

    int translateToSensorRefSP(const int spIndexRefID);
    int translateToIndexRefSP(const int spSensorRefID);

    void displayBinBoundaries();
    void displayBinWeightDistribution();
    void displayBinWeightDistributionLayer(int layer);
    void closePlot();

    void testbinSearch();
    void testCoordImpl();
    void testSPID();
};


#endif //COSMICTRIGGER_PATTERNENGINE_H
