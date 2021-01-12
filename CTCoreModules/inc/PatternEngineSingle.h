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
#include "SPCalculations.h"


class PatternEngineSingle : public SPCalculations{
    /**
     * Pattern Engine Single handles one single station of the detector and its super pixel mapping.
     * This includes the creation of cartesian decision borders between super pixels, Super Pixel Frequency monitoring,
     * Cartesian hit -> Super Pixel assignment.
     *
     * It was build with several mapping systematics in mind, that is the purpose of the ``mode`` parameter, to switch
     * between these.
     * Forseen was a radial SP mapping and a planar x-projective mapping.
     * Completely implemented and used was only the former. To use a projective mappaing, some functions have to be
     * implemented still.
     *
     * The reason for the SP borders being stored in an array is the idea of a non-uniform mapping where not all super
     * pixels are evenly sized. This could be a means to increase bkg discrimination power in areas with a lot of background.
     * The Pattern Engine could be upgraded to a non-uniform mapping, just by adding a memberfunction of a separate
     * constructor, that sets up the bin-boundary vecotors in the desired way.
     *
     * Furhter Comment: If the plotting functions are used, remember to call closePlot(); before deleting the
     * Pattern Engine, as otherwise the pdf output wont be readable.
     */

private:
    int getLayer(const float r);
    int getZSP(const float z);
    int getPhiSP(const float phi);
    int getXSP(const float x);
    int getSPZcoord(const int sp2DID);
    int getSPWcoord(const int sp2DID);

    std::string getAreaTag();

    template <typename T>
    int binSearch(std::vector<T> targetVector, T value, int start, int end);

    template <typename T>
    int seqSearch(std::vector<T> targetVector, T value);

    void displayBinWeightDistributionLayer(int layer);

    float getRadius(const float x, const float y);
    float getHitPhi(const float x, const float y);

    void initializeSPbins();
    void initializeMembers(const int spWpartition, const int spZpartition, const int area, const int mode);

    std::string pltf();

    //necessary for wbin mapping with projective x-coord, not phi
    float layerFactor[4] = {0.0}; // scale down x boundaries with reference to outer layer

    int mode = 0; //default is radial sp (phi as variable): mode = 0
    int area = 0;
    int wBinCount{}; //width (phi or x)
    int zBinCount{};
    int totalBinCount{};

    //for the region the PES manages
    float centralDetectorZmin;
    float centralDetectorZmax;
    float centralDetectorLength;

    float widthMax{};
    float widthMin{};
    float widthDimLength{};

    std::vector<float> phiBins;
    std::vector<float> zBins;
    std::vector<float> xBins; //not used so far -> x-projective Mapping

    float SPZsize{};
    float SPWsize{};

    //count freq. of hits in each super pixel
    std::vector<int> spWeights[4];

    std::string plottingpath = "plots/PatternEngineSingle/";
    std::string plottingfile = "PEplots.pdf"; //be careful! reset in constructor with area
    std::string runspecs;
    int firstplot = 1;

public:
    PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode /*, const int area*/);
    PatternEngineSingle(const int spWpartition, const int spZpartition, const int mode, /*const int area,*/ std::string plottingpath);
    PatternEngineSingle(const int spWpartition, const int spZpartition, const int mode, const int area, std::string plottingpath);

    PatternEngineSingle();
    //~PatternEngineSingle(); //Improve: Move closePlot() to destructor.

    unsigned int getSuperPixel(const float x, const float y, const float z);

    // not implemented yet, but important for later implementation on FPGA -->
    int translateToSensorRefSP(const int spIndexRefID);
    int translateToIndexRefSP(const int spSensorRefID);
    // <--

    void displayBinBoundaries();
    void displayBinWeightDistribution();
    std::string getRunSpecs();
    void closePlot(); //IMPORTANT! plots are all written to one pdf, which is closed by this function.

    //tests
    void testbinSearch();
    void testCoordImpl();
    void testSPID();
};


#endif //COSMICTRIGGER_PATTERNENGINESINGLE_H
