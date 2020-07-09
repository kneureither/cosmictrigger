//
// Created by Konstantin Neureither on 16.06.20.
//

#include <vector>
#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <assert.h>
#include <iostream>

#include "utilityFunctions.h"
#include "plots.h"
#include "PatternEngineSingle.h"
#include "SuperPixelHit.h"
#include "basicDefines.h"

PatternEngineSingle::PatternEngineSingle() = default;

//initialize with default path
PatternEngineSingle::PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode) {
    printf("----INIT SINGLE PATTERN ENGINE----\n");
    this->initializeMembers(spXpartition, spZpartition, 0, mode);
    this->initializeSPbins();
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PE INITIATED---\n\n");

}

//initalize as single PE
PatternEngineSingle::PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode, std::string plottingpath) {
    printf("----INIT SINGLE PATTERN ENGINE----\n");
    this->initializeMembers(spXpartition, spZpartition, 0, mode);
    this->initializeSPbins();
    this->plottingpath = plottingpath;
    this->runspecs = "wbins" + get_string(spXpartition) + "zbins" + get_string(spZpartition);
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PE INITIATED---\n\n");
}

//initialize as one detector segment of compound PE
PatternEngineSingle::PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode, const int area, std::string plottingpath) {
    printf("----INIT SINGLE PATTERN ENGINE AS PARTIAL BY PATTERN ENGINE----\n");
    printf("Zmin=%f Zmax=%f\n", DetectorZmin[area], DetectorZmax[area]);
    this->centralDetectorZmin = DetectorZmin[area];
    this->centralDetectorZmax = DetectorZmax[area];
    this->initializeMembers(spXpartition, spZpartition, area, mode);
    this->initializeSPbins();
    this->plottingpath = plottingpath;
    this->plottingfile = "PEplots_" + getAreaTag() + ".pdf";
    this->runspecs = "wbins" + get_string(spXpartition) + "zbins" + get_string(spZpartition);
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PES INITIATED---\n\n");
}

//init members (called by constructors)
void PatternEngineSingle::initializeMembers(const int spXpartition, const int spZpartition, const int area, const int mode) {

    this->wBinCount = spXpartition;
    this->zBinCount = spZpartition;
    this->totalBinCount = spXpartition * spZpartition;
    this->centralDetectorLength = this->centralDetectorZmax - this->centralDetectorZmin;
    this->widthMin = -PI;
    this->widthMax = PI;
    this->widthDimLength = this->widthMax - this->widthMin;
    this->SPZsize = centralDetectorLength / zBinCount;
    this->SPWsize = widthDimLength / wBinCount;
    this->mode = mode;
    this->area = area;

    printf("z SP size=%f phi SP size=%f\n", SPZsize, SPWsize);
    printf("wBinCount=%d zbinCount=%d totalBinCount=%d\n", wBinCount, zBinCount, totalBinCount);
    printf("wdithMin=%f widthMax=%f widthDimLength=%f\n", widthMin, widthMax, widthDimLength);

    for(int i=0; i<4; i++) {
//        this->spWeights[i].reserve(this->totalBinCount);
        for(int j=0; j<(this->totalBinCount); j++) this->spWeights[i].push_back(0);
        printf("initialized spWeights[%d]\n", i);
    }
}

//init the superpixel bin boundary arrays
void PatternEngineSingle::initializeSPbins() {
    float SPZlength = centralDetectorLength / zBinCount;
    for(int i=0; i <= zBinCount; i++) {
        zBins.push_back(i * SPZlength + centralDetectorZmin);
    }

    float SPWlength = widthDimLength / wBinCount;
    for(int i = 0; i <= wBinCount; i++) {
        phiBins.push_back(i * SPWlength + widthMin);
    }

    //TODO implement xBins
}

////SOME UTILITY FUNCTIONS ------------------->
int PatternEngineSingle::getZSP(const float z) {
    return binSearch(this->zBins, z, 0, zBins.size() - 1);
//    return seqSearch(this->zBins, z);

}

int PatternEngineSingle::getPhiSP(const float phi) {
    return binSearch(this->phiBins, phi, 0, this->phiBins.size() - 1);
//    return seqSearch(this->phiBins, phi);
}


int PatternEngineSingle::getLayer(const float r) {
    if (layerBoundaries[0] <= r && r <= layerBoundaries[1]) {
        return 0;
    } else if (layerBoundaries[1] < r && r <= layerBoundaries[2]) {
        return 1;
    } else if (layerBoundaries[2] < r && r <= layerBoundaries[3]) {
        return 2;
    } else if (layerBoundaries[3] < r && r <= layerBoundaries[4]) {
        return 3;
    } else {
        return -1;
    }
}

float PatternEngineSingle::getRadius(const float x, const float y) {
    return sqrt(pow(x,2) + pow(y,2));
}

float PatternEngineSingle::getHitPhi(const float x, const float y) {
    return atan2(y, x);
}

int PatternEngineSingle::getSPXcoord(const int sp2DID) {
    return sp2DID % this->wBinCount;
}

int PatternEngineSingle::getSPZcoord(const int sp2DID) {
    return sp2DID / this->wBinCount;
}

template<typename T>
int PatternEngineSingle::binSearch(std::vector<T> targetVector, T value, int start, int end) {
    int size = end-start;
    if (size <= 0) {
        return -1;
    }
    if(size == 1) {
        return start;
    }
    int center = (start + end) / 2;
    if (targetVector[center] <= value && value <= targetVector[end]) {
        return binSearch(targetVector, value, center, end);
    } else if (targetVector[start] <= value && value < targetVector[center]){
        return binSearch(targetVector, value, start, center);
    } else {
        return -1;
    }
}

template<typename T>
int PatternEngineSingle::seqSearch(std::vector<T> targetVector, T value) {
    for(int i = 0; i < targetVector.size() - 1; i++) {
        if(targetVector[i] <= value && value < targetVector[i+1]){
            return i;
        }
    }
    if(value == targetVector[targetVector.size()-1]) {
        return targetVector.size()-2;
    } else {
        return -1;
    }
}

//// <-------------------- SOME UTILITY FUNCTIONS


//get xyz coordinates and calculate the superpixel from it. most important function
unsigned int PatternEngineSingle::getSuperPixel(const float x, const float y, const float z) {
//    printf("Coordinates delivered: x=%f, y=%f, z=%f\n", x,y,z);
    float r = getRadius(x, y);
//    printf("radius=%f\n", r);
    int layer = getLayer(r);
//    printf("layer=%d\n", layer);
    int phiSPIndex = getPhiSP(getHitPhi(x, y));
//    printf("phiSP=%d\n", phiSPIndex);
    int zSPIndex = getZSP(z);
//    printf("zSP=%d\n", zSPIndex);

    assert(0 <= area && area < 3);
    assert(0 <= layer && layer < 4);
    assert(totalBinCount < 4096); // fatal because SID needs format of 0xFFF

    int sp2DID = computeIndex(zSPIndex, phiSPIndex, this->wBinCount);
    if(sp2DID < 0 || totalBinCount < sp2DID) sp2DID = 0;
    unsigned int SPID = computeSPID(layer, area, sp2DID);

//    int SPID = layer + 10 * area + 100 * sp2DID;
    printf("Hit SuperPixel params\t x=%f, y=%f, z=%f, layer=%d zIndex=%d, phiIndex=%d, Sp2D=%d, SPID=%#X\n", x,y,z,layer, zSPIndex, phiSPIndex, sp2DID, SPID);

    //    printf("SP2D = %d, weights = %d \n", sp2DID, this->spWeights[0][sp2DID]);
    if(SPID != 0) this->spWeights[layer].at(sp2DID) = this->spWeights[layer].at(sp2DID) + 1;

    return SPID;
}

//// TESTS ( ONLY WORK ON CENTRAL DETECTOR SO FAR) FIXME ------------------->
void PatternEngineSingle::testbinSearch() {
    std::vector<float> bins;
    printf("Testing binSearch().... \n Added data: \n");
    for(int i = 0; i < 20; i++) {
        bins.push_back(i*PI);
        printf("  %f", i*PI);
    }
    int result;
    float data;

    printf("\n\nIn-Between test data results \n");
    for(int i = 0; i<19; i++) {
        data = PI/2 + i*PI;
        result = binSearch(bins, data, 0, bins.size() - 1);
        printf("Testing data=%f, index result=%d\n", data, result);
        assert(i==result);
    }

    printf("\n\nSpot on test data results \n");
    for(int i = 0; i<19; i++) {
        data = i*PI;
        result = binSearch(bins, data, 0, bins.size() - 1);
        printf("Testing data=%f, index result=%d\n", data, result);

        assert(i==result);
    }

    printf("\n\nOut of bounds test \n");
    data = -PI;
    result = binSearch(bins, data, 0, bins.size() - 1);
    printf("Testing data=%f, index result=%d\n", data, result);
    assert(result == -1);

    data = 21 * PI;
    result = binSearch(bins, data, 0, bins.size() - 1);
    printf("Testing data=%f, index result=%d\n", data, result);
    assert(result == -1);
}

void PatternEngineSingle::testCoordImpl() {
    int z = -100;
    int y = 10;
    int x = 10;

    float phi = getHitPhi(x, y);
    int phiSP = getPhiSP(phi);
    int zSP = getZSP(z);

    int SP2D = zSP * wBinCount + phiSP;
    assert(getLayer(getRadius(x,y)) == 0);
    assert(phiSP == getSPXcoord(SP2D));
    assert(zSP == getSPZcoord(SP2D));
}

void PatternEngineSingle::testSPID() {
    unsigned int sp1 = getSuperPixel(1, 0, -100);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));

    sp1 = getSuperPixel(0, 0, -200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(20, 0, 0);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(45, 0, 100);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(55, 0, 200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(65, 0, 200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));

    float x=55.0;
    float y=3.0;
    float z=100;

    unsigned int sp2=getSuperPixel(x,y,z);
    int phiSP = getPhiSP(getHitPhi(x, y));
    int zSP = getZSP(z);
    int index = computeIndex(zSP, phiSP, this->wBinCount);
    int layer = getLayer(getRadius(x,y));
    int area=0;
    printf("-- custom evaluation: x=%f, y=%f, z=%f, phiSP=%d, zSP=%d, SPindex=%d, layer=%d, area=%d\n",x,y,z,phiSP, zSP, index, layer, area);
    unsigned int SPID = computeSPID(layer, area, index);
    int reslayer = getLayerFromSPID(SPID);
    int resarea = getAreaFromSPID(SPID);
    int resindex = getIndexFromSPID(SPID);
    printf("-- result for SPID dec=%u hex=%#X    ->   layer=%d, area=%d, index=%d\n", SPID, SPID, reslayer, resarea, resindex);

    assert(reslayer == layer);
    assert(resindex == index);
    assert(resarea == area);

//    getIndexFromSPID(0x1234F000); //assertion must fail
}
//// <------------  TESTS ( ONLY WORK ON CENTRAL DETECTOR SO FAR)



//// PRINT GRAPHS AND PLOTS ----------->
void PatternEngineSingle::displayBinBoundaries() {

    auto *canvas = new TCanvas("bin_boundaries", "bin_boundaries", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * grid = new TH2F("h", ("Superpixel bin grid boundaries "+ this->getAreaTag()).c_str() , zBinCount, centralDetectorZmin, centralDetectorZmax, wBinCount, widthMin, widthMax);
    grid->SetStats(false);
    grid->GetYaxis()->SetBinLabel(1, "-#pi");
    grid->GetYaxis()->SetBinLabel((int) (grid->GetNbinsY() / 2), "0");
    grid->GetYaxis()->SetBinLabel(grid->GetNbinsY(), "#pi");
    labelAxis(grid, "z [mm]", "#phi", 0.07);
    grid->Draw();

    std::vector<TLine*>boundaries;

    for(int i=0; i < phiBins.size(); i++) {
        TLine * l = new TLine(this->centralDetectorZmin, this->phiBins[i], this->centralDetectorZmax, this->phiBins[i]);
        l->SetLineColor(kBlue);
        boundaries.push_back(l);
    }

    for(int i=0; i < zBins.size(); i++) {
        TLine * l = new TLine(this->zBins[i], widthMin, this->zBins[i], widthMax);
        l->SetLineColor(kRed);
        boundaries.push_back(l);
    }

    for(int i = 0; i< boundaries.size(); i++) boundaries[i]->Draw();

    canvas->Print((this->plottingpath + this->runspecs + pltf()).c_str(), "pdf");

}

void PatternEngineSingle::displayBinWeightDistributionLayer(int layer) {
    std::string layerno = get_string(layer);
    auto *canvas = new TCanvas(("bin_weights"+layerno).c_str(), ("bin_weights"+layerno).c_str(), 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * h_binweights = new TH2F("h", ("Superpixel bin weights layer " + layerno + " " + getAreaTag()).c_str(), zBinCount, 0, zBinCount, wBinCount, 0, wBinCount);
    h_binweights->SetStats(true);
    labelAxis(h_binweights, "z bins", "#phi bins");

    //Fill 2d hist with weights
    for(int i=1; i<totalBinCount; i++) {
        for(int j=1; j<spWeights[layer].at(i); j++) {
            h_binweights->Fill(getSPZcoord(i+1), getSPXcoord(i+1));
        }
    }

    // Norm
    h_binweights -> Scale(1.0 / h_binweights->Integral());

//    std::cout << "Bin z :" << h_binweights->GetXaxis()->GetBinWidth(1) << std::endl;

    h_binweights->Draw("colz");

    canvas->Print((this->plottingpath + this->runspecs + pltf()).c_str(), "pdf");

    TFile * tF = new TFile("output.root", "RECREATE");
    h_binweights->Write();
    TH1F * tP = (TH1F*)h_binweights->ProjectionY();
    tP->Write();
    tF->Close();

}

void PatternEngineSingle::displayBinWeightDistribution() {
    for(int i=0; i<4; i++) {
        this->displayBinWeightDistributionLayer(i);
    }
}

std::string PatternEngineSingle::pltf() {
    if(this->firstplot == 1) {
        printf("first file\n");
        this->firstplot = 0;
        return plottingfile + "(";
    } else {
        return plottingfile;
    }
}

void PatternEngineSingle::closePlot() {
    auto *c_final = new TCanvas("c_final", "c_final");
    c_final->Print((this->plottingpath + this->runspecs + pltf() + ")").c_str(), "pdf");
}

std::string PatternEngineSingle::getAreaTag() {
    if(this->area == 0) {
        return "central";
    } else if (this->area == 1) {
        return "recurlR";
    } else if (this->area == 2) {
        return "recurlL";
    } else {
        return "area" + get_string(this->area);
    }
}
//// <------------ PRINT GRAPHS AND PLOTS





