//
// Created by Konstantin Neureither on 16.06.20.
//
#ifndef PI
#define PI 3.1415926535
#endif //PI

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
#include "PatternEngine.h"
#include "../util/plots.h"
#include "SuperPixelHit.h"
#include "../util/utilityFunctions.h"

PatternEngine::PatternEngine(float spXpartition, float spZpartition, int mode) {
    printf("----INIT PATTERN ENGINE----\n");
    this->initializeMembers(spXpartition, spZpartition, mode);
    this->initializeSPbins();
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PE INITIATED---\n\n");

}

PatternEngine::PatternEngine(float spXpartition, float spZpartition, int mode, std::string plottingpath) {
    printf("----INIT PATTERN ENGINE----\n");
    this->initializeMembers(spXpartition, spZpartition, mode);
    this->initializeSPbins();
    this->plottingpath = plottingpath;
    this->runspecs = "wbins" + get_string(spXpartition) + "zbins" + get_string(spZpartition);
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PE INITIATED---\n\n");
}

void PatternEngine::initializeMembers(float spXpartition, float spZpartition, int mode) {

    this->wBinCount = spXpartition;
    this->zBinCount = spZpartition;
    this->totalBinCount = spXpartition * spZpartition;
    this->centralDetectorLength = this->centralDetectorZmax - this->centralDetectorZmin;
    this->widthMin = -PI;
    this->widthMax = PI;
    this->widthDimLength = this->widthMax - this->widthMin;
    this->SPZsize = centralDetectorLength / zBinCount;
    this->SPPhisize = widthDimLength / wBinCount;

    printf("z SP size=%f phi SP size=%f\n", SPZsize, SPPhisize);
    printf("wBinCount=%d zbinCount=%d totalBinCount=%d\n", wBinCount, zBinCount, totalBinCount);
    printf("wdithMin=%f widthMax=%f widthDimLength=%f\n", widthMin, widthMax, widthDimLength);

    for(int i=0; i<4; i++) {
//        this->spWeights[i].reserve(this->totalBinCount);
        for(int j=0; j<(this->totalBinCount); j++) this->spWeights[i].push_back(0);
        printf("initialized spWeights[%d]\n", i);
    }
}

void PatternEngine::initializeSPbins() {
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

int PatternEngine::getXSP(float x) {
    return 0;
}

int PatternEngine::getZSP(float z) {
    return binSearch(this->zBins, z, 0, zBins.size() - 1);
//    return seqSearch(this->zBins, z);

}

int PatternEngine::getPhiSP(float phi) {
    return binSearch(this->phiBins, phi, 0, this->phiBins.size() - 1);
//    return seqSearch(this->phiBins, phi);
}


int PatternEngine::getLayer(float r) {
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

int PatternEngine::getSuperPixel(float x, float y, float z) {
//    printf("Coordinates delivered: x=%f, y=%f, z=%f\n", x,y,z);
    float r = getRadius(x, y);
//    printf("radius=%f\n", r);
    int layer = getLayer(r);
//    printf("layer=%d\n", layer);
    int phiSPIndex = getPhiSP(getHitPhi(x, y));
//    printf("phiSP=%d\n", phiSPIndex);
    int zSPIndex = getZSP(z);
//    printf("zSP=%d\n", zSPIndex);

    int area = 0;

    int sp2DID = zSPIndex * wBinCount + phiSPIndex;

    if(sp2DID < 0 || totalBinCount < sp2DID) sp2DID = 0;

    int SPID = layer + 10 * area + 100 * sp2DID;
    //printf("Hit SuperPixel params\t x=%f, y=%f, z=%f, layer=%d zIndex=%d, phiIndex=%d, Sp2D=%d, SPID=%d\n", x,y,z,layer, zSPIndex, phiSPIndex, sp2DID, SPID);

    //    printf("SP2D = %d, weights = %d \n", sp2DID, this->spWeights[0][sp2DID]);
    if(SPID != 0) this->spWeights[layer].at(sp2DID) = this->spWeights[layer].at(sp2DID) + 1;

    return SPID;
}

float PatternEngine::getRadius(float x, float y) {
    return sqrt(pow(x,2) + pow(y,2));
}

float PatternEngine::getHitPhi(float x, float y) {
    return atan2(y, x);
}

int PatternEngine::getSPXcoord(int sp2DID) {
    return sp2DID % this->wBinCount;
}

int PatternEngine::getSPZcoord(int sp2DID) {
    return sp2DID / this->wBinCount;
}

template<typename T>
int PatternEngine::binSearch(std::vector<T> targetVector, T value, int start, int end) {
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
int PatternEngine::seqSearch(std::vector<T> targetVector, T value) {
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

void PatternEngine::testbinSearch() {
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

void PatternEngine::testCoordImpl() {
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


void PatternEngine::displayBinBoundaries() {

    auto *canvas = new TCanvas("bin_boundaries", "bin_boundaries", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * grid = new TH2F("h", "Superpixel bin grid boundaries", zBinCount, centralDetectorZmin, centralDetectorZmax, wBinCount, widthMin, widthMax);
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

void PatternEngine::displayBinWeightDistributionLayer(int layer) {
    std::string layerno = get_string(layer);
    auto *canvas = new TCanvas(("bin_weights"+layerno).c_str(), ("bin_weights"+layerno).c_str(), 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * h_binweights = new TH2F("h", ("Superpixel bin weights layer " + layerno).c_str(), zBinCount, 0, zBinCount, wBinCount, 0, wBinCount);
    h_binweights->SetStats(true);
    labelAxis(h_binweights, "z bins", "#phi bins");


    for(int i=1; i<totalBinCount; i++) {
        for(int j=1; j<spWeights[layer].at(i); j++) {
            h_binweights->Fill(getSPZcoord(i+1), getSPXcoord(i+1));
        }
    }

    // Norm
    h_binweights -> Scale(1.0 / h_binweights->Integral());


    std::cout << "Bin z :" << h_binweights->GetXaxis()->GetBinWidth(1) << std::endl;


    h_binweights->Draw("colz");

    canvas->Print((this->plottingpath + this->runspecs + pltf()).c_str(), "pdf");

    TFile * tF = new TFile("output.root", "RECREATE");
    h_binweights->Write();
    TH1F * tP = (TH1F*)h_binweights->ProjectionY();
    tP->Write();
    tF->Close();

}

void PatternEngine::displayBinWeightDistribution() {
    for(int i=0; i<4; i++) {
        this->displayBinWeightDistributionLayer(i);
    }
}

std::string PatternEngine::pltf() {
    if(this->firstplot == 1) {
        printf("first file\n");
        this->firstplot = 0;
        return plottingfile + "(";
    } else {
        return plottingfile;
    }
}

void PatternEngine::closePlot() {
    auto *c_final = new TCanvas("c_final", "c_final");
    c_final->Print((this->plottingpath + this->runspecs + pltf() + ")").c_str(), "pdf");

}





