//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_PATTERNENGINESINGLE_H
#define COSMICTRIGGER_PATTERNENGINESINGLE_H

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
#include "SPCalculations.h"


//TODO Implement other mapping of superpixels (maybe just a cartesian one
// -- This could either be done by deriving cartesian and phi dependent from a common base class
// -- Or include it into this class and use the int mode parameter


class PatternEngineSingle : public SPCalculations{

private:
    int getLayer(const float r);
    int getZSP(const float z);
    int getPhiSP(const float phi);
    int getSPZcoord(const int sp2DID);
    int getSPXcoord(const int sp2DID);

    std::string getAreaTag();

    template <typename T>
    int binSearch(std::vector<T> targetVector, T value, int start, int end);

    template <typename T>
    int seqSearch(std::vector<T> targetVector, T value);

    void displayBinWeightDistributionLayer(int layer);

    float getRadius(const float x, const float y);
    float getHitPhi(const float x, const float y);

    void initializeSPbins();
    void initializeMembers(const int spXpartition, const int spZpartition, const int area, const int mode);

    std::string pltf();


    float layerFactor[4] = {0.0}; // scale down x boundaries with reference to outer layer

    int mode = 0; //default is cylindrical mode = 0
    int area = 0;
    int wBinCount{}; //width (phi or x)
    int zBinCount{};
    int totalBinCount{};

    //only center region implemented for now
    float centralDetectorZmin;
    float centralDetectorZmax;
    float centralDetectorLength;

    float widthMax{};
    float widthMin{};
    float widthDimLength{};

    std::vector<float> phiBins;
    std::vector<float> zBins;
    std::vector<float> xBins;

    float SPZsize{};
    float SPWsize{};

    std::vector<int> spWeights[4];

    std::string plottingpath = "plots/PatternEngineSingle/";
    std::string plottingfile = "PEplots.pdf"; //be careful! reset in constructor with area
    std::string runspecs;
    int firstplot = 1;

public:
    PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode /*, const int area*/);
    PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode, /*const int area,*/ std::string plottingpath);
    PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode, const int area, std::string plottingpath);

    PatternEngineSingle();
//    ~PatternEngineSingle();

    unsigned int getSuperPixel(const float x, const float y, const float z);

    int translateToSensorRefSP(const int spIndexRefID);
    int translateToIndexRefSP(const int spSensorRefID);

    void displayBinBoundaries();
    void displayBinWeightDistribution();
    std::string getRunSpecs();
    void closePlot();

    void testbinSearch();
    void testCoordImpl();
    void testSPID();

};


#endif //COSMICTRIGGER_PATTERNENGINESINGLE_H
